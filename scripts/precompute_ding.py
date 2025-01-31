import data_utils as du
from argparse import ArgumentParser
import random as r
import os
import subprocess as sp
import sys


DING_HOME=os.path.join(os.path.dirname(os.path.dirname(sys.argv[0])),"tools","ding-cf-main")


GUROBI_CMD = 'gurobi_cl'

def pairs_for_mode(mode,n,leaves):
    pairs = []
    print("Precomputing for",len(leaves),"leaves")
    if mode=='linear':
        print("Mode: linear")
        leaves_left=list(leaves)
        while len(leaves_left)>1:
            l1=r.choice(leaves_left)
            leaves_left.remove(l1)
            l2= r.choice(leaves_left)
            leaves_left.remove(l2)
            pairs.append((l1,l2))
        if len(leaves_left)==1:
            l1=leaves_left.pop()
            l2 =r.choice([l for l in leaves if l != l1])
        return pairs
    pairs = []
    for l1 in leaves:
        for l2 in leaves:
            if l1 < l2:
                pairs.append((l1,l2))
    if mode=='all':
        print("Mode: all")
        return pairs
    else:
        print("Mode: random")
        return list(r.sample(pairs,max(n,len(pairs))))
    

def cp_tree(edges):
    tree = {}
    for child, parent in edges:
        assert(child not in tree)
        tree[child]=parent
    return tree

def lca_trace_cp_tree(tree,a,b):
    a_trace_ = [a]
    curr = a
    while curr in tree:
        curr = tree[curr]
        a_trace_.append(curr)
    
    b_trace = [b]
    curr = b
    a_seen = set(a_trace_)
    while curr not in a_seen:
        curr=tree[curr]
        b_trace.append(curr)
    a_trace = a_trace_[:a_trace_.index(curr)]
    print(a_trace,b_trace,curr)
    return b_trace+a_trace

MODES = ['all','linear','random']

parser = ArgumentParser()
parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
parser.add_argument('workdir')
parser.add_argument('tree', type=open,
            help='phylogenetic tree as parent-child relation table')
parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
parser.add_argument('--family-bounds','-fmb',type=open)

parser.add_argument('--mode',choices=MODES,default='all')
parser.add_argument('--total-timelimit',type=int,default=10*60)
parser.add_argument('-n',type=int,default=1)
parser.add_argument("--ding",help="Path to the ding home directory",default=DING_HOME)
parser.add_argument("--gurobi",default=GUROBI_CMD)
parser.add_argument("outfile")
args = parser.parse_args()



run_ding=os.path.join(args.ding,'ding_cf.py')

candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
    sep=args.separator)

families = candidateAdjacencies['families']
speciesTree = du.parseTree(args.tree)
leaves = set([x for x, v in du.getLeaves(speciesTree).items() if v])

pairs = pairs_for_mode(args.mode,args.n,leaves)

print("Number pairs",len(pairs))

tree = cp_tree(speciesTree)




print(tree)

if args.family_bounds:
    fam_bounds = du.parseFamilyBounds(args.family_bounds)
else:
    fam_bounds = dict()
du.fillFamilyBounds(candidateAdjacencies['families'],fam_bounds)

unimogs = dict()
for l in leaves:
    unimogs[l]=du.get_unimog(l,candidateAdjacencies['genes'],candidateAdjacencies['adjacencies'],args.separator)


best_bounds = {}
timelim = int(args.total_timelimit/len(pairs))
for a,b in pairs:
    fam_max = dict()
    fam_min = dict()
    file_prefix = "{a}_{b}".format(a=a,b=b)
    genomes = lca_trace_cp_tree(tree,a,b)
    for g in genomes:
        for f in fam_bounds[g]:
            if f not in fam_max:
                fam_max[f]=fam_bounds[g][f][1]
                fam_min[f]=fam_bounds[g][f][0]
            fam_min[f]=min(fam_min[f],fam_bounds[g][f][0])
            #use minimum here because we will have to sort via a genome with
            #this number of markers
            fam_max[f]=min(fam_max[f],fam_bounds[g][f][1])
    modelfilename = os.path.join(args.workdir,file_prefix+".mmodel")
    unimogfilename = os.path.join(args.workdir,file_prefix+'.unimog')
    #write dict to model file
    with open(modelfilename,'w') as mmfile:
        for f in fam_max:
            print("\t".join([f,str(fam_min[f]),str(fam_max[f])]),file=mmfile)
    #write unimog file
    with open(unimogfilename,'w') as ugfile:
        print(unimogs[a],file=ugfile)
        print(unimogs[b],file=ugfile)
    #call ding
    ilpfilename=os.path.join(args.workdir,file_prefix+'.ilp')
    ilplog = os.path.join(args.workdir,file_prefix+'.ilp.log')
    with open(ilplog,'w') as log:
        sp.run(["python3",run_ding,unimogfilename,"-c",modelfilename,"--writeilp",ilpfilename],stderr=log)
    #call gurobi
    gurobilog = os.path.join(args.workdir,file_prefix+'.gurobi.log')
    with open(gurobilog,"w") as log:
        sp.run([args.gurobi,"Threads=1","TimeLimit={}".format(timelim),ilpfilename],stdout=log)
    #parse gurobi logfile
    with open(gurobilog) as log:
        lastline = log.readlines()[-1]
        bbs = [bb.strip() for bb in lastline.split(",") if bb.strip().startswith("best bound")]
        if len(bbs)!=1:
            print("Corrupted best bound file, skipping...")
        else:
            bbstr=bbs.pop()
            bb=float(bbstr[len("best bound "):])
            best_bounds[(a,b)]=bb
#output into distance list
with open(args.outfile,"w") as f:
    for (a,b),bb in best_bounds.items():
        print(a,b,bb,file=f)


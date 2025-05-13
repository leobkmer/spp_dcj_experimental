#!/usr/bin/env python3

import data_utils as du
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('tree', type=open,
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
parser.add_argument('family_bounds',help="Marker ranges per a")
parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')

args = parser.parse_args()

with open(args.family_bounds) as fmb:
    fam_bounds=du.parseFamilyBounds(fmb)

candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)
speciesTree = du.parseTree(args.tree)

tree = du.cp_tree(speciesTree)


def weighted_adj_freq_in_subtree(tree,candidateAdjacencies,fbounds,leaves,sep):
    fam_adj_freqs = {}
    kill_adjacencies = {}
    for l in leaves:
        #print(candidateAdjacencies['adjacencies'])
        fam_adj_freqs[l]={}
        kill_adjacencies[l]=set()
        for x in candidateAdjacencies['adjacencies'][l]:
            
            (a,x),(b,y)=x
            a=a.split(sep)[0]
            b=b.split(sep)[0]
            adj = (a,x),(b,y)
            adj = tuple(sorted(adj))
            fam_adj_freqs[l][adj]=1
            if len(candidateAdjacencies['families'][l][a])==1 and len(candidateAdjacencies['families'][l][b])==1:
                kill_adjacencies[l].add(adj)
    root,pc_tree = du.cp_to_pc(tree)
    traversal = [root]
    curr_level = [root]
    while True:
        next_level = []
        for x in curr_level:
            next_level.extend(pc_tree.get(x,[]))
        traversal.extend(next_level)
        curr_level=next_level
        if len(curr_level)==0:
            break
    #bottom up traversal
    traversal=traversal[::-1]
    
    for v in traversal:
        if not v in pc_tree:
            continue
        fam_adj_freqs[v]={}
        kill_adjacencies[v] = set.intersection(*[set(((a,x),(b,y)) for ((a,x),(b,y)) in  kill_adjacencies[child] 
                                                   if fbounds[v].get(a,None)==(1,1) and fbounds[v].get(b,None)==(1,1))
                                                   for child in pc_tree[v]])
        for child in pc_tree[v]:
            assert(child in fam_adj_freqs)
            for xtr,frq in fam_adj_freqs[child].items():
                if not xtr in fam_adj_freqs[v]:
                    fam_adj_freqs[v][xtr]=0
                fam_adj_freqs[v][xtr]+=frq/len(pc_tree[v])
    return fam_adj_freqs,kill_adjacencies

#def all_possible_circular_adjacencies(gene_dict):
#    all_genes = set()
#    for genes in gene_dict.values():
#        for g in genes:
#            all_genes.add(g)
#    all_extremities = set([(g,du.EXTR_HEAD) for g in all_genes])
#    all_extremities=all_extremities.union(set([(g,du.EXTR_TAIL) for g in all_genes]))
#    all_adjacencies = [(x,y) for x in all_extremities for y in all_extremities if x < y]
#    return all_adjacencies


weights = candidateAdjacencies['weights']
leaves=set([x for x, v in du.getLeaves(speciesTree).items() if v])

freqs,kill_adjacencies = weighted_adj_freq_in_subtree(tree,candidateAdjacencies,fam_bounds,leaves,args.separator)
#for x in (freqs.values()):
#    print(x.items())

#adj_inner = all_possible_circular_adjacencies(candidateAdjacencies['genes'])




def all_adj_from_bound(fbounds,kill_adjacencies,separator):
    poss_anc = dict()
    
    for genome, fams in fbounds.items():
        all_genes = set()
        kill_extremities = set()
        rem_adj = []
        for ((a,x),(b,y)) in kill_adjacencies[genome]:
            e1=("{}{}1".format(a,separator),x)
            e2=("{}{}1".format(b,separator),y)
            e1,e2 = (e1,e2) if e1 < e2 else (e2,e1)
            rem_adj.append((e1,e2))
            kill_extremities.add(e1)
            kill_extremities.add(e2)
        print("Kill adjacencies for  {}: {}".format(genome,kill_adjacencies[genome]),file=sys.stderr)
        for fname,(_,high) in fams.items():
            all_genes.update(["{fname}{sep}{i}".format(fname=fname,i=i,sep=separator) for i in range(1,high+1)])
        all_extremities = set([(g,du.EXTR_HEAD) for g in all_genes])
        all_extremities.update([(g,du.EXTR_TAIL) for g in all_genes])
    
        all_adjacencies = [(x,y) for x in all_extremities for y in all_extremities if x < y and (x not in kill_extremities) and (y not in kill_extremities)] 
        poss_anc[genome] = all_adjacencies+rem_adj
    return poss_anc



root,pc_tree = du.cp_to_pc(tree)

possible_adjacencies = all_adj_from_bound(fam_bounds,kill_adjacencies,separator=args.separator)

#add adjacencies for inner nodes
for anc in pc_tree:
    for x_,y_ in possible_adjacencies[anc]:
        (a,x),(b,y)=x_,y_
        a=a.split(args.separator)[0]
        b=b.split(args.separator)[0]
        fadj = (a,x),(b,y)
        fadj = tuple(sorted(fadj))
        if not fadj in freqs[anc]:
            weights[anc][(x_,y_)]=0
        else:
            weights[anc][(x_,y_)]=freqs[anc][fadj]

print("#Species\tGene_1\tExt_1\tSpecies\tGene_2\tExt_2\tWeight")
for gnm in weights:
    for ((a,x),(b,y)),w in weights[gnm].items():
        print("\t".join([gnm,a,x,gnm,b,y,str(w)]))
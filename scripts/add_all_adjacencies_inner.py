#!/usr/bin/env python3

import data_utils as du
from argparse import ArgumentParser
import sys



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



def weighted_adj_freq(tree,adjacencies,fbounds,leaves,sep):
    fam_adj_freqs = {}
    #kill_adjacencies are selected according to Lemma 7
    kill_adjacencies = {}
    #semi_kill_adjacencies are selected according to Lemma 8
    semi_kill_adjacencies = {}
    for l in leaves:
        fam_adj_freqs[l]={}
        kill_adjacencies[l]=set()
        semi_kill_adjacencies[l]=set()
        leaf_extremity_counts=dict()
        for x in adjacencies[l]:
            a,b=x
            if not a in leaf_extremity_counts:
                leaf_extremity_counts[a]=0
            if not b in leaf_extremity_counts:
                leaf_extremity_counts[b]=0
            leaf_extremity_counts[a]+=1
            leaf_extremity_counts[b]+=1
        for x in adjacencies[l]:
            (a,x),(b,y)=x
            isaunique = leaf_extremity_counts[(a,x)]==1
            isbunique = leaf_extremity_counts[(b,y)]==1
            a=a.split(sep)[0]
            b=b.split(sep)[0]
            adj = (a,x),(b,y)
            adj = tuple(sorted(adj))
            fam_adj_freqs[l][adj]=1
            if fbounds[l].get(a,None)==(1,1) and fbounds[l].get(b,None)==(1,1) and isaunique and isbunique:
                kill_adjacencies[l].add(adj)
    root,pc_tree = du.cp_to_pc(tree)
    traversal = bottom_up_traversal(root, pc_tree)
    for v in traversal:
        if not v in pc_tree:
            continue
        fam_adj_freqs[v]={}
        kill_adjacencies[v] = set.intersection(*[set(((a,x),(b,y)) for ((a,x),(b,y)) in  kill_adjacencies[child] 
                                                   if fbounds[v].get(a,None)==(1,1) and fbounds[v].get(b,None)==(1,1))
                                                   for child in pc_tree[v]])
        semi_kill_adjacencies[v]=set()
        if len(pc_tree[v])==2:
            #Lemma 8 works only for exactly two children
            pot_semi = dict()
            for child in pc_tree[v]:
                pot_semi[child]=dict()
                for ((a,x),(b,y)) in  kill_adjacencies[child]:                             
                    if fbounds[v].get(a,None)==(1,1) and fbounds[v].get(b,None)==(1,1):
                        if ((a,x),(b,y)) not in kill_adjacencies[v]:
                            pot_semi[child][(a,x)]=(b,y)
                            pot_semi[child][(b,y)]=(a,x)
            c1 = pc_tree[v][0]
            c2 = pc_tree[v][1]
            handled = set()
            for x in pot_semi[c1]:
                if x in handled:
                    continue
                y = pot_semi[c1][x]
                if not y in pot_semi[c2]:
                    continue
                z = pot_semi[c2][y]
                if not x in pot_semi[c2]:
                    continue
                w = pot_semi[c2][x]
                if x in handled or y in handled or z in handled or w in handled:
                    continue
                if pot_semi[c1].get(z,None)==w:
                    #4 cycle found, Lemma 8 applies
                    semi_kill_adjacencies[v].update([(x,y),(w,z),(w,x),(y,z)])
                handled.update([x,y,z,w])
        for child in pc_tree[v]:
            assert(child in fam_adj_freqs)
            for xtr,frq in fam_adj_freqs[child].items():
                if not xtr in fam_adj_freqs[v]:
                    fam_adj_freqs[v][xtr]=0
                fam_adj_freqs[v][xtr]+=frq/len(pc_tree[v])
    return fam_adj_freqs,kill_adjacencies#,semi_kill_adjacencies

def bottom_up_traversal(root, pc_tree):
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
    return traversal
    

#freqs,kill_adjacencies = weighted_adj_freq_in_subtree(tree,candidateAdjacencies,fam_bounds,leaves,args.separator)
#for x in (freqs.values()):
#    print(x.items())

#adj_inner = all_possible_circular_adjacencies(candidateAdjacencies['genes'])




def all_adj_from_bound(fbounds,kill_adjacencies,separator,prove_filter=True):
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
        #for ((a,x),(b,y)) in semi_kill_adjacencies[genome]:
        #    e1=("{}{}1".format(a,separator),x)
        #    e2=("{}{}1".format(b,separator),y)
        #    e1,e2 = (e1,e2) if e1 < e2 else (e2,e1)
        #    rem_adj.append((e1,e2))
        #    kill_extremities.add(e1)
        #    kill_extremities.add(e2)
        #print("Kill adjacencies for  {}: {}".format(genome,kill_adjacencies[genome]),file=sys.stderr)
        for fname,(_,high) in fams.items():
            all_genes.update(["{fname}{sep}{i}".format(fname=fname,i=i,sep=separator) for i in range(1,high+1)])
        all_extremities = set([(g,du.EXTR_HEAD) for g in all_genes])
        all_extremities.update([(g,du.EXTR_TAIL) for g in all_genes])
    
        all_adjacencies = [(x,y) for x in all_extremities for y in all_extremities if x < y and ((x not in kill_extremities) and (y not in kill_extremities) or not prove_filter)] 
        poss_anc[genome] = list(set(all_adjacencies+rem_adj))
    return poss_anc




def invert_mrk(mrk):
    r = []
    for s,m in mrk:
        r.append(("+" if s=="-" else "-",m))
    return r[::-1]


def canonicalize_chrm(chrm):
    lc,mrk = chrm
    if lc == du.CHR_LINEAR:
        return (du.CHR_LINEAR,min(mrk,invert_mrk(mrk)))
    assert(lc==du.CHR_CIRCULAR)
    #circular chromosome
    #inefficient, but for simplicity's sake
    rotations=[mrk[i:]+mrk[:i] for i in range(len(mrk))]
    rrotatiions = [invert_mrk(y) for y in rotations]
    return (du.CHR_CIRCULAR,min(rotations+rrotatiions))

def canonicalize_gnm(gnm):
    name,chrs = gnm
    can_chrs = [canonicalize_chrm(chrm) for chrm in chrs]
    return (name,sorted(can_chrs,key=lambda c: (sorted([m for _,m in c[1]]),c)))





def canon_adjacencies(gnm):
    return [tuple(sorted(adj)) for adj in du.unimog2adjacencies(canonicalize_gnm(gnm))]


def canonicalized_has_distance0or1(adja,adjb):
    #check if we find a distance of 0 or 1 based on canonical adjacency representations
    #Note: not all distance 1 cases can be found this way
    diffa = set(adja).difference(set(adjb))
    diffb = set(adjb).difference(set(adja))
    #print(diffa)
    #print(diffb)
    if len(diffa)==0 and len(diffb)==0:
        return 0
    if len(diffa)==1 and len(diffb)==2:
        axs=list(sorted(diffa.pop()))
        bxs = list(sorted([(m,x) for y in diffb for (m,x) in y if x!= du.EXTR_CAP]))
        if axs==bxs:
            return 1
    if len(diffb)==1 and len(diffa)==2:
        bxs=list(sorted(diffb.pop()))
        axs = list(sorted([(m,x) for y in diffa for (m,x) in y if x!= du.EXTR_CAP]))
        if axs==bxs:
            return 1
        

    if len(diffa)!=2:
        return None
    if len(diffb)!=2:
        return None
    ((a,b),(c,d))=tuple(sorted([sorted(x) for x in diffa]))
    ((e,f),(g,h))=tuple(sorted([sorted(x) for x in diffb]))
    if not a==e:
        return None
    if c==f and tuple(sorted((b,d)))==(g,h):
        return 1
    if d==f and tuple(sorted((b,c)))==(g,h):
        return 1
    return None


def add_unimog_to_marker_ranges(fam_bounds,unimog):
    for name, chrs in unimog:
        assert(name not in fam_bounds)
        fam_bounds[name]=dict()
        for _,mrk in chrs:
            for _,m in mrk:
                (x,y)=fam_bounds[name].get(m,(0,0))
                assert(x==y)
                fam_bounds[name][m]=(x+1,y+1)



UNIVERSE = "all"
def trim_tree_fitch(tree,leaves,unimog):
    filtered_nodes = []
    root,pc_tree=du.cp_to_pc(tree)
    possible_linearizations = dict()
    for gnm in unimog:
        name,_=gnm
        ca = canon_adjacencies(gnm)
        possible_linearizations[name]=[ca]
    traversal = bottom_up_traversal(root,pc_tree)
    for name in traversal:
        if name in leaves:
            continue
        children = pc_tree[name]
        if len(children) > 2:
            possible_linearizations[name]=UNIVERSE
        if len(children)==1:
            possible_linearizations[name]=possible_linearizations[children[0]]
        if len(children)==2:
            c1 = children[0]
            c2 = children[1]
            if possible_linearizations[c1]==UNIVERSE or possible_linearizations[c2]==UNIVERSE:
                possible_linearizations[name]=UNIVERSE
            else:
                fitchable = True
                has_intersection = False
                for a in possible_linearizations[c1]:
                    for b in possible_linearizations[c2]:
                        d = canonicalized_has_distance0or1(a,b)
                        if d is None:
                            fitchable=False
                        if d==0:
                            has_intersection=True
                if not fitchable:
                    possible_linearizations[name]=UNIVERSE
                elif has_intersection:
                    filtered_nodes.append(c1)
                    filtered_nodes.append(c2)
                    possible_linearizations[name]=[x for x in possible_linearizations[c1] if x in possible_linearizations[c2]]
                else:
                    filtered_nodes.append(c1)
                    filtered_nodes.append(c2)
                    possible_linearizations[name]=possible_linearizations[c1]+possible_linearizations[c2]
    return filtered_nodes,possible_linearizations


parser = ArgumentParser()
parser.add_argument('tree',
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
parser.add_argument('unimog',
            help='unimog file for the leaf genomes')
parser.add_argument('family_bounds',help="Marker ranges in the tree")
parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
parser.add_argument('--no-proof-filter',action='store_true')
parser.add_argument('--no-tree-trimming',action='store_true')
parser.add_argument("--write-adjacencies")
parser.add_argument("--write-tree")
args = parser.parse_args()

with open(args.family_bounds) as fmb:
    fam_bounds=du.parseFamilyBounds(fmb)


with open(args.unimog) as ug:
    unimog = du.parseUniMoG(ug)




canon_ad1 = canon_adjacencies(unimog[1])
canon_ad2 = canon_adjacencies(unimog[2])


add_unimog_to_marker_ranges(fam_bounds,unimog)


#adjacencies = dict()
#for genome in unimog:
#    name,_ = genome
#    adj = du.unimog2adjacencies(genome)
#    adjacencies[name]=adj


with open(args.tree) as trf:
    speciesTree = du.parseTree(trf)


tree = du.cp_tree(speciesTree)
root,cptr = du.cp_to_pc(tree)
rootchild = cptr[root][0]

leaves=set([x for x, v in du.getLeaves(speciesTree).items() if v])


if not args.no_tree_trimming:
    fltn,lins=trim_tree_fitch(tree,leaves,unimog)
    fltn=set(fltn)
    print("Able to filter {} nodes from the tree.".format(len(fltn)),file=sys.stderr)

    for f in fltn:
        if f==rootchild:
            continue
        del tree[f]
        del fam_bounds[f]
    
    if len(tree)>1 and rootchild in fltn:
        del tree[rootchild]
        del fam_bounds[rootchild]
    leaves=set([x for x, v in du.getLeaves(list(tree.items())).items() if v])

    print("New leaves",leaves,file=sys.stderr)
    adjacencies = dict()
    for gname, listofadjsets in lins.items():
        if gname not in tree:
            continue
        if listofadjsets==UNIVERSE:
            continue
        adjacencies[gname]=set()
        for adj in listofadjsets:
            adjacencies[gname].update(adj)
else:
    adjacencies = dict()
    for genome in unimog:
        name,_ = genome
        adj = du.unimog2adjacencies(genome)
        adjacencies[name]=adj 


print("Propagating adjacency weights through tree..",file=sys.stderr)


freqs, kills = weighted_adj_freq(tree,adjacencies,fam_bounds,leaves,sep=args.separator)

#print(semikills,file=sys.stderr)

inner_bounds = dict()
for gnm,bnds in fam_bounds.items():
    if gnm not in leaves:
        inner_bounds[gnm]=bnds

root, pct = du.cp_to_pc(tree)


if len(pct[root])==2:
    #root is not necessary, remove and connect the children directly
    c1,c2 = tuple(pct[root])
    del inner_bounds[root]
    del kills[root]
    #del semikills[root]
    tree[c2]=c1
    del tree[c1]
    print("Removed unneccesary root",root,"new tree: ",tree,file=sys.stderr)


adjs = all_adj_from_bound(inner_bounds,kills,args.separator,prove_filter=not args.no_proof_filter)

for gnm, adj in adjs.items():
    assert(gnm not in adjacencies)
    adjacencies[gnm]=adj

if args.write_adjacencies:
    with open(args.write_adjacencies,"w") as adjf:
        for gnm,adj in adjacencies.items():
            ws = freqs.get(gnm,dict())
            for (a,x),(b,y) in adj:
                af=a.split(args.separator)[0]
                bf=b.split(args.separator)[0]
                adj = tuple(sorted(((af,x),(bf,y))))
                w = ws.get(adj,0)
                print(f"{gnm}\t{a}\t{x}\t{gnm}\t{b}\t{y}\t{w}",file=adjf)

if args.write_tree:
    with open(args.write_tree,"w") as tf:
        for child, parent in tree.items():
            print(f"{child}\t{parent}",file=tf)
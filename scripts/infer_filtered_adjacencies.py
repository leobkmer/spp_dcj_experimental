from tree_trim import trim_tree_fitch, UNIVERSE
from argparse import ArgumentParser
import data_utils as du







def find_filtered_subtrees(pc_tree,root,filtered):
    active_nodes = [root]
    filtered_roots = []
    while len(active_nodes)>0:
        n = active_nodes.pop()
        assert(n not in filtered)
        filter_children = set()
        unfilter_children = set()
        if n not in pc_tree:
            continue
        for child in pc_tree[n]:
            if child in filtered:
                filter_children.add(child)
            else:
                unfilter_children.add(child)
        assert(len(filter_children)==0 or len(unfilter_children)==0)
        if len(filter_children)==0:
            for child in unfilter_children:
                active_nodes.append(child)
        else:
            filtered_roots.append(n)
    return filtered_roots



def set_filter_adjacencies(cptr,root,root_adj,linearizations):
    set_adjacencies = dict()
    set_adjacencies[root]=set(tuple(sorted(a)) for a in root_adj)
    active_nodes = [root]
    while len(active_nodes)>0:
        n = active_nodes.pop()
        n_adj = set_adjacencies[n]
        for child in cptr[n]:
            preferred_lin = None
            for lin in linearizations[child]:
                #print(child,lin)
                lin = set(lin)
                preferred_lin=lin
                if lin == n_adj:
                    break
            set_adjacencies[child]=preferred_lin
            if child in cptr:
                active_nodes.append(child)
    return set_adjacencies


    
    





parser = ArgumentParser()
parser.add_argument('tree',
            help='original phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
parser.add_argument('unimog',
            help='unimog file for the leaf genomes')
parser.add_argument("adj",help="Resolved adjacencies file")

args = parser.parse_args()

with open(args.unimog) as ug:
    unimog = du.parseUniMoG(ug)


with open(args.adj) as f:
    adjacencies = du.parseAdjacencies(f.readlines())["adjacencies"]
    




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

fltn,lins=trim_tree_fitch(tree,leaves,unimog)

for froot in find_filtered_subtrees(cptr,root,fltn):
    for name, adjs in set_filter_adjacencies(cptr,froot,adjacencies[froot],lins).items():
        if name==froot or name==rootchild:
            continue
        for (a,x),(b,y) in adjs:
            print(f"{name}\t{a}\t{x}\t{name}\t{b}\t{y}\t1")
#!/usr/bin/env python3

import data_utils as du
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('tree', type=open,
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')

args = parser.parse_args()


candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)
speciesTree = du.parseTree(args.tree)

tree = du.cp_tree(speciesTree)


def infer_fam_adj_freqs(tree,candidateAdjacencies,leaves,sep):
    fam_adj_freqs = {}
    for l in leaves:
        #print(candidateAdjacencies['adjacencies'])
        fam_adj_freqs[l]={}
        for x in candidateAdjacencies['adjacencies'][l]:
            
            (a,x),(b,y)=x
            a=a.split(sep)[0]
            b=b.split(sep)[0]
            adj = (a,x),(b,y)
            adj = tuple(sorted(adj))
            fam_adj_freqs[l][adj]=1
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
        for child in pc_tree[v]:
            assert(child in fam_adj_freqs)
            for xtr,frq in fam_adj_freqs[l].items():
                if not xtr in fam_adj_freqs[v]:
                    fam_adj_freqs[v][xtr]=0
                fam_adj_freqs[v][xtr]+=frq/len(pc_tree[v])
    return fam_adj_freqs

def all_possible_circular_adjacencies(gene_dict):
    all_genes = set()
    for genes in gene_dict.values():
        for g in genes:
            all_genes.add(g)
    all_extremities = set([(g,du.EXTR_HEAD) for g in all_genes])
    all_extremities=all_extremities.union(set([(g,du.EXTR_TAIL) for g in all_genes]))
    all_adjacencies = [(x,y) for x in all_extremities for y in all_extremities if x < y]
    return all_adjacencies


weights = candidateAdjacencies['weights']
leaves=set([x for x, v in du.getLeaves(speciesTree).items() if v])

freqs = infer_fam_adj_freqs(tree,candidateAdjacencies,leaves,args.separator)

adj_inner = all_possible_circular_adjacencies(candidateAdjacencies['genes'])


root,pc_tree = du.cp_to_pc(tree)
#add adjacencies for inner nodes
for anc in pc_tree:
    for x_,y_ in adj_inner:
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
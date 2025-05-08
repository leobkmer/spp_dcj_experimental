#!/usr/bin/env python3


from argparse import ArgumentParser
from Bio import Phylo

parser = ArgumentParser()
parser.add_argument("nwk")
parser.add_argument("out")

args = parser.parse_args()


tree = Phylo.read(args.nwk, "newick")
i=0
for clade in tree.find_clades():
    if clade.name is None:
        clade.name="anc%d"%i
        i+=1


Phylo.write(tree,args.out,"newick")
import csv
import sys
from argparse import ArgumentParser
import data_utils as du

#genome gene_disambiguator group
def read_groups(gf):
    gm = dict()
    with open(gf) as f:
        for entries in csv.reader(f,delimiter='\t'):
            genome,geneuq,grp = entries
            if not grp in gm:
                gm[grp]=[]
            gm[grp].append((genome,geneuq))
    return gm


def cleanup_groups(gm):
    gm_clean = dict()
    gid = 0
    for gn,g in gm:
        gid+=1
        has_paralogs=False
        gnm_grp = dict()
        for gnm,geneuq in g:
            if not gnm in gnm_grp:
                gnm_grp[gnm] = []
            else:
                has_paralogs = True
            gnm_grp[gnm].append(geneuq)
        if has_paralogs:
            print("Warning: group {} contains paralogs. It will be arbitrarily split".format(gn),file=sys.stderr)
            gnm_remaining = set(gnm_grp.values())
            while len(gnm_remaining) > 0:
                gnm_exiting = set()
                gid+=1
                gm_clean[gid]=[]
                for x in gnm_remaining:
                    front = gnm_grp[x].pop()
                    gm_clean[gid].append((x,front))
                    gm_clean[gid].append()
                    if len(gnm_grp[x])==0:
                        gnm_exiting.add(x)
        else:
            gm_clean[gid]=g
    return gid


def adjacency_fitch(tree,adjacencies):
    pass



def main():
    parser = ArgumentParser()
    parser.add_argument('tree',
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
    parser.add_argument('candidateAdjacencies',
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('mgroups',help='File of pre-grouped markers. If groups are not fully resolved, they will be arbitrarily split.')
    args = parser.parse_args()
    with open(args.candidateAdjacencies) as cA:
        candAdj = du.parseAdjacencies(cA,sep=args.separator)
    with open(args.tree) as treefile:
        speciesTree = du.parseTree(treefile)
    mgroups = cleanup_groups(read_groups(args.mgroups))


main()

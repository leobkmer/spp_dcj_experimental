import data_utils as du
from argparse import ArgumentParser

parser = ArgumentParser()


parser.add_argument('tree', type=open,
            help='phylogenetic tree as parent-child relation table')
parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
args = parser.parse_args()


speciesTree = du.parseTree(args.tree)
leaves = set([x for x, v in du.getLeaves(speciesTree).items() if v])
candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)
bounds = dict()
fams = candidateAdjacencies['families']
for gnm in fams:
    if gnm in leaves:
        continue
    bounds[gnm]=dict()
    for f in fams[gnm]:
        bounds[gnm][f]=(1,len(fams[gnm][f]))

for gnm in bounds:
    for f in bounds[gnm]:
        print("\t".join([gnm,f,str(bounds[gnm][f][0]),str(bounds[gnm][f][1])]))
import os
from argparse import ArgumentParser


def read_genomes(gnmdir):
    genomes = []
    for fl in os.listdir(gnmdir):
        print(fl)
        if not fl.endswith("_GENOME.tsv"):
            continue
        print("reading file",fl)
        name = fl.split("_GENOME.tsv")[0]
        filename = os.path.join(gnmdir,fl)
        with open(filename) as f:
            genes_raw = [line.split() for line in f.readlines()[1::]]
            genes = []
            for g in genes_raw:
                g[0] = int(g[0])
                print(g)
                genes.append(tuple(g))
            genes.sort()
            genes_no_pos=[(o,fam,subid) for _,fam,o,subid in genes]
            genomes.append((name,genes_no_pos))
    return genomes


def genomes2adjacencies(genomes):
    adj = []
    for name, genes in genomes:
        o,fam,subid = genes[-1]

        last = fam+"_"+subid
        laste = "h" if o=="+" else "t"
        for (o,fam,subid) in genes:
            thise = "t" if o=="+" else "h"
            this = fam+"_"+subid
            adj.append((name,last,laste,name,this,thise,"0"))
            last = this
            laste = "h" if o=="+" else "t"
    print(adj)
    return adj


def main():
    parser = ArgumentParser()
    parser.add_argument('genomeDir')
    parser.add_argument("outfile")
    args = parser.parse_args()
    genomes = read_genomes(args.genomeDir)
    with open(args.outfile,"w") as f:
        for a in genomes2adjacencies(genomes):
            print("\t".join(a),file=f)
    

main()
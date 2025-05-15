#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
import csv
import re

# import from own packages
import data_utils as du
from collections import defaultdict


def get_f_values(solfilen):
    fvs = {}
    with open(solfilen) as solf:
        for line in solf:
            if line.startswith("#"):
                continue
            n,v= line.split()
            if n.startswith("f") and "_" in n:
                fvs[n.split('_')[1]]=float(v)
    return fvs  


def get_name_map(peif):
    nm = {}
    with open(peif) as pf:
        for line in pf:
            a,b, n = line.split()
            nm[n]=(a,b)
    return nm

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sol_file',
            help='solution file of GUROBI optimizer')
    parser.add_argument('id_to_extremity_map',type=open,
            help='mapping between node IDs and extremities')
    parser.add_argument("--edge-distances")
    parser.add_argument("--phylogeny-edge-ids")
    args = parser.parse_args()

    #
    # load data
    #

    # id-to-extremity mapping
    id2extr,_ = du.parse_id_file(args.id_to_extremity_map)#dict()
    id2ext = dict([(str(v),k) for k,v in id2extr.items()])
    #for line in csv.reader(args.id_to_extremity_map, delimiter = '\t'):
    #    if line:
    #        id2ext[line[0]] = tuple(line[1:])
    with open(args.sol_file) as sf:
        adjacenciesList,_ = du.parseSOLAdj(sf, id2ext)
    # write adjacencies
    du.writeAdjacencies(adjacenciesList,defaultdict(lambda: defaultdict(float)), stdout)
    if args.edge_distances:
        with open(args.edge_distances,"w") as f:
                if args.phylogeny_edge_ids:
                        nm = get_name_map(args.phylogeny_edge_ids)
                        to_id = lambda x: nm[x]
                else:
                        to_id = lambda x : x
                for eid, v in get_f_values(args.sol_file).items():
                        a,b = to_id(eid)
                        print(a,b,v,file=f)

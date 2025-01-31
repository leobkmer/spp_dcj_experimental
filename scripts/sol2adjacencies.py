#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
import csv
import re

# import from own packages
import data_utils as du
from collections import defaultdict

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sol_file', type=open,
            help='solution file of GUROBI optimizer')
    parser.add_argument('id_to_extremity_map', type=open,
            help='mapping between node IDs and extremities')

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

    adjacenciesList,_ = du.parseSOLAdj(args.sol_file, id2ext)
    # write adjacencies
    du.writeAdjacencies(adjacenciesList,defaultdict(lambda: defaultdict(float)), stdout)

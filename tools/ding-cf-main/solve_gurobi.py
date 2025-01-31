#!/usr/bin/env python3.7
import gurobipy as gp
from gurobipy import GRB
import sys
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("lpfile", help="File containing the (I)LP")
parser.add_argument("solfile", help="File to write solution to.")
parser.add_argument("-t",default=1,type=int, help="Number of threads to use.")
parser.add_argument("--timelim",type=int, default=60*60*24, help="Timelimit in seconds.")
args = parser.parse_args()
model = gp.read(args.lpfile)
model.Params.Threads = args.t
model.Params.Timelimit = args.timelim
model.optimize()
model.write(args.solfile)

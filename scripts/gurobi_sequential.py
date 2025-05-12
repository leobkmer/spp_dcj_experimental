#!/usr/bin/env python3.7
import gurobipy as gp
import sys
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("lpfile", help="File containing the (I)LP")
parser.add_argument("solfile", help="File to write solution to.")
parser.add_argument("-t",default=1,type=int, help="Number of threads to use.")
parser.add_argument("--tree",help="File containing the tree")
parser.add_argument("--warm-start")
#parser.add_argument("--timelim",type=int, default=60*60*24, help="Timelimit in seconds.")
args = parser.parse_args()
model = gp.read(args.lpfile)
model.Params.Threads = args.t
w=model.getVarByName("w")
f=model.getVarByName("f")
if args.warm_start:
    model.read(args.warm_start)
    model.update()
for alph in [0.1,0.9]:
    model.Params.Timelimit = 500
    ialph = alph -1
    model.setObjective(ialph*w+alph*f)
    model.update()
    print("Optimizing with objective: ",model.getObjective())
    model.optimize()
model.setObjective(f)
model.Params.Timelimit = 2600
model.update()
print("------------------FINAL OPTIMIZATION----------------")
model.optimize()
model.write(args.solfile)
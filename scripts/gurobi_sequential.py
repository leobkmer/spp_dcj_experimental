import gurobipy as gp
import sys
from argparse import ArgumentParser
import data_utils as du





parser = ArgumentParser()

parser.add_argument("lpfile", help="File containing the (I)LP")
parser.add_argument("solfile", help="File to write solution to.")
parser.add_argument("-t",default=1,type=int, help="Number of threads to use.")
parser.add_argument("--weightsolves",nargs="*",type=float,default=[0.1,0.5,0.9],help="Solve a weighted version of the ILP first.")
parser.add_argument("--weightsolve-proportion",type=float,default=0.2,help="Maximum proportion of time spent on the weighted versions.")
parser.add_argument("--timelim",type=int,default=3600,help="Overall time limit in seconds.")
parser.add_argument("--memlim",type=int,default=10000,help="Overall mem limit in MB.")
parser.add_argument("--tree",help="File containing the tree")
parser.add_argument("--warm-start")
#parser.add_argument("--timelim",type=int, default=60*60*24, help="Timelimit in seconds.")
args = parser.parse_args()
if args.tree:
    tree_ids = du.read_tree_edge_name_map(args.tree)
    root, pc_tree = du.cp_to_pc(du.cp_tree(tree_ids.keys()))
    

model = gp.read(args.lpfile)
model.Params.Threads = args.t
model.Params.SoftMemLimit = args.memlim/1000 - 5
model.Params.Timelimit = args.timelim
w=model.getVarByName("w")
f=model.getVarByName("f")
if args.warm_start:
    model.read(args.warm_start)
    model.update()
model.Params.MIPFocus=1
model.Params.MIPGap=1/100
for alph in args.weightsolves:
    tl = args.weightsolve_proportion/len(args.weightsolves)*args.timelim
    model.Params.Timelimit = tl
    ialph = alph -1
    model.setObjective(ialph*w+alph*f)
    model.update()
    print("-------------------WEIGHTED PRE-OPTIMIZATION:Optimizing with objective: ",model.getObjective()," and time limit {} ----------------------".format(tl))
    model.optimize()
    #print("Intermediate solution with f-objective: {}".format(f.X))
model.setObjective(f)
model.Params.MIPFocus=0
model.Params.MIPGap=1e-4
model.update()
print("------------------FINAL OPTIMIZATION----------------")
model.optimize()
model.write(args.solfile)
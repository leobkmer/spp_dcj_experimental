import gurobipy as gp
import sys
from argparse import ArgumentParser
import data_utils as du
import random as r
from math import floor



parser = ArgumentParser()

parser.add_argument("lpfile", help="File containing the (I)LP")
parser.add_argument("solfile", help="File to write solution to.")
parser.add_argument("tree",help="File containing the tree")
parser.add_argument("-t",default=1,type=int, help="Number of threads to use.")
parser.add_argument("--timelim",type=int,default=3600,help="Overall time limit in seconds.")
parser.add_argument("--weightsolves",nargs="*",type=float,default=[0.1,0.5,0.9],help="Solve a weighted version of the ILP first.")
parser.add_argument("--weightsolve-proportion",type=float,default=0.3,help="Maximum proportion of time spent on the weighted versions.")
parser.add_argument("--treesolve-proportion",type=float,default=0.3,help="Maximum proportion of time spent on solving subtrees first.")
parser.add_argument("--min-subtree-size",type=int,default=3)
parser.add_argument("--max-subtree-size",type=int,default=25)
parser.add_argument("--subsample-subtrees",help="Subsample subtrees to this number",type=int)
parser.add_argument("--warm-start")
args = parser.parse_args()
model = gp.read(args.lpfile)
model.Params.Threads = args.t
w=model.getVarByName("w")
f=model.getVarByName("f")
if args.warm_start:
    model.read(args.warm_start)
    model.update()





tree_ids_r = du.read_tree_edge_name_map(args.tree)


model.Params.MIPFocus=1
model.Params.MIPGap=1/100
for alph in [0]+args.weightsolves:
    tl = max(args.weightsolve_proportion/(len(args.weightsolves)+1)*args.timelim,100)
    model.Params.Timelimit = tl
    ialph = alph -1
    model.setObjective(ialph*w+alph*f)
    model.update()
    print("-------------------WEIGHTED PRE-OPTIMIZATION with objective: ",model.getObjective()," and time limit {} ----------------------".format(tl))
    model.optimize()
    
    if alph==0:
        max_effect_of_weights=abs(model.ObjBound)

for x in model.getVars():
    try:
        model.getVarByName(x.VarName).setAttr("Start",x.X)
    except AttributeError:
        break
tree_ids = {}
for k,v in tree_ids_r.items():
    tree_ids[v]=k
root, pc_tree = du.cp_to_pc(du.cp_tree(tree_ids.keys()))
subtrees = [du.get_subtree_edges(x,pc_tree) for x in du.subtrees_ascending_size(root,pc_tree) if x!=root]
subtrees = [s for s in subtrees if len(s)>= args.min_subtree_size and len(s) <= args.max_subtree_size]
if args.subsample_subtrees and len(subtrees) > args.subsample_subtrees:
    subtrees=sorted(list(r.choices(subtrees,k=args.subsample_subtrees)),key=lambda x: len(x))
model.setObjective(f)

if len(subtrees) > 0:
    st_timelim=max(100,args.timelim*args.treesolve_proportion/(len(subtrees)))
else:
    st_timelim=0
print("Optimizing lower bound by solving ",len(subtrees), "subtrees")
model.Params.MIPFocus=2
model.update()
for i,subtree in enumerate(subtrees,start=1):    
    tree_edge_ids = [tree_ids[e] for e in subtree]
    fvars = [model.getVarByName("f_{}".format(i)) for i in tree_edge_ids]
    wvar = model.getVarByName("w")
    model.setObjective(sum(fvars))#-1/(max_effect_of_weights+1)*wvar)
    model.Params.Timelimit = st_timelim
    model.update()
    print("---------------------------TREE PRE-OPTIMIZATION for subtree {}/{} ({}) ".format(i,len(subtrees),subtree),"with objective: ",model.getObjective()," and time limit {} ----------------------".format(st_timelim))
    model.optimize()
    model.addConstr(sum([model.getVarByName("f_{}".format(i)) for i in tree_edge_ids])>=floor(model.ObjBound))
    for x in model.getVars():
        if not "_" in x.VarName:
            continue
        v_id_split = x.VarName.split("_")
        if  v_id_split[0].startswith("r") or  v_id_split[0].startswith("p") or v_id_split[0] in ["x","f","w","l","z","y","p","d","w","b","n","c","s"]:
            if v_id_split[1] in tree_edge_ids:
                try:
                    model.getVarByName(x.VarName).setAttr("VarHintVal",x.X)
                except AttributeError:
                    break
    model.update()




model.setObjective(f)
model.Params.MIPFocus=0
model.Params.MIPGap=1/100
model.update()
full_timelim=max(100,(1-args.weightsolve_proportion - args.treesolve_proportion)*args.timelim)
model.Params.Timelimit = args.timelim

print("------------------FINAL OPTIMIZATION----------------")
model.optimize()
model.write(args.solfile)
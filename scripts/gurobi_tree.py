import gurobipy as gp
import sys
from argparse import ArgumentParser
import data_utils as du





parser = ArgumentParser()

parser.add_argument("lpfile", help="File containing the (I)LP")
parser.add_argument("solfile", help="File to write solution to.")
parser.add_argument("tree",help="File containing the tree")
parser.add_argument("-t",default=1,type=int, help="Number of threads to use.")
parser.add_argument("--timelim",type=int,default=3600,help="Overall time limit in seconds.")
parser.add_argument("--memlim",type=int,default=10000,help="Overall mem limit in MB.")
args = parser.parse_args()
model = gp.read(args.lpfile)
model.Params.Threads = args.t
model.Params.SoftMemLimit = args.memlim/1000 - 5
w=model.getVarByName("w")
f=model.getVarByName("f")
#if args.warm_start:
#    model.read(args.warm_start)
#    model.update()


tree_ids_r = du.read_tree_edge_name_map(args.tree)
tree_ids = {}
for k,v in tree_ids_r.items():
    tree_ids[v]=k
root, pc_tree = du.cp_to_pc(du.cp_tree(tree_ids.keys()))
subtrees = [du.get_subtree_edges(x,pc_tree) for x in du.subdivide_tree(root,pc_tree,max_leaves=5)]
subtrees = [subtree for subtree in subtrees if len(subtree) > 1]

for i,x in enumerate(subtrees):
    for y in subtrees[i+1::]:
        assert(len(set(x).intersection(set(y)))==0)
        t_idsx = set([tree_ids[a] for a in x])
        t_idsy = set([tree_ids[a] for a in y])
        assert(len(t_idsx.intersection(t_idsy))==0)
        v_idsx = set([a for a,_ in x]+[b for _,b in x])
        v_idsy = set([a for a,_ in y]+[b for _,b in y])
        assert(len(v_idsx.intersection(v_idsy))==0)

#remainder = set(du.get_subtree_edges(root,pc_tree))


st_timelim=max(100,args.timelim/(len(subtrees))*0.5)
print("Optimizing first over ",len(subtrees), "subtrees")
remainder_model = model.copy()
for i,subtree in enumerate(subtrees,start=1):
    print("---------------------------Optimizing for subtree {} ({}/{})---------------------------".format(subtree,i,len(subtrees)))
    m=model.copy()
    tree_edge_ids = [tree_ids[e] for e in subtree]
    fvars = [m.getVarByName("f_{}".format(i)) for i in tree_edge_ids]
    m.setObjective(sum(fvars))
    m.Params.Timelimit = st_timelim
    m.update()
    m.optimize()
    model.addConstr(sum([model.getVarByName("f_{}".format(i)) for i in tree_edge_ids])>=m.ObjBound)
    for x in m.getVars():
        if not "_" in x.VarName:
            continue
        v_id_split = x.VarName.split("_")
        if  v_id_split[0].startswith("r") or  v_id_split[0].startswith("p") or v_id_split[0] in ["x","f","w","l","z","y","p","d","w","b","n","c","s"]:
            if v_id_split[1] in tree_edge_ids:
                try:
                    model.getVarByName(x.VarName).setAttr("Start",x.X)
                except AttributeError:
                    break

    model.update()


#remainder_model.update()
#remainder_model.setObjective(remainder_model.getVarByName("f"))
#remainder_model.Params.Timelimit = 200
#print("----------------------Optimizing remainder-------------------")
#remainder_model.update()
#remainder_model.optimize()

#for x in remainder_model.getVars():
#       model.getVarByName(x.VarName).setAttr("Start",x.X)



model.setObjective(f)
model.Params.Timelimit = args.timelim
#model.Params.MIPFocus=0
#model.Params.MIPGap=1e-4
model.update()
print("------------------FINAL OPTIMIZATION----------------")
model.optimize()
model.write(args.solfile)
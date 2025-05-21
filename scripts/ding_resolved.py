from argparse import ArgumentParser
import os
import subprocess as sp

parser = ArgumentParser()
parser.add_argument("unimog")
parser.add_argument("tree")
parser.add_argument("workdir")

args = parser.parse_args()

DING_PATH =  '../../tools/ding-cf-main/ding_cf.py'
GUROBI = 'gurobi_cl'
edges = []
with open(args.tree) as tf:
    tf.readline()
    for line in tf:
        child,parent = line.strip().split()
        edges.append((child,parent))

branch_lengths=dict()
for child, parent in edges:
    ilpfile = os.path.join(args.workdir,child+"_"+parent+".ilp")
    solfile = os.path.join(args.workdir,child+"_"+parent+".sol")
    sp.run(["python3",DING_PATH,args.unimog,"-p",child,parent,"--writeilp",ilpfile])
    sp.run([GUROBI,"ResultFile={}".format(solfile),ilpfile])
    with open(solfile) as sf:
        for line in sf:
            if not "Objective value" in line:
                    continue
            sv = float(line.strip().split("=")[1])
        branch_lengths[(child,parent)]=sv

finalfile = os.path.join(args.workdir,"result.txt")
tot_len = sum(branch_lengths.values())
with open(finalfile,"w") as f:
    print(tot_len,file=f)
    for c,p in sorted(branch_lengths.keys()):
        print(c,p,branch_lengths[(c,p)],file=f)
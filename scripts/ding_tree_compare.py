from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument("ding")
parser.add_argument("spp2")
args = parser.parse_args()

dists = []
with open(args.ding) as f:
    dv = float(f.readline().strip())
    for line in f:
        dists.append(float(line.strip().split()[-1]))

fvars = []
with open(args.spp2) as s2f:
    for line in s2f:
                if line.startswith("f"):
                    if not "_" in line:
                         continue
                    fvar, x= line.strip().split()
                    f,n = fvar.split("_")
                    fvars.append(((f,n),float(x)))
                if not "Objective value" in line:
                    continue
                s2v = float(line.strip().split("=")[1])

print(dv,s2v)
ok=True

if abs(dv-s2v) >= 0.001:
     ok=False

for u,fv in zip(dists,sorted(fvars)): 
    f,v = fv
    print(f,u,v)
    if abs(v-u) >= 0.001:
         ok=False

if not ok:
    sys.exit(1)
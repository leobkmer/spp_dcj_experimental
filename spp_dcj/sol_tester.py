#!/usr/bin/env python3

from argparse import ArgumentParser
import sys
parser = ArgumentParser()
parser.add_argument("solfile")
parser.add_argument("--solval",type=float)
parser.add_argument("--vars",type=str,nargs='*',default=[])
parser.add_argument("--vals",type=int,nargs='*',default=[])

args = parser.parse_args()



valdict = dict(zip(args.vars,args.vals))

ok =True
with open(args.solfile) as f:
    obj_txt = '# Objective value = '
    for line in f:
        if line.startswith(obj_txt):
            obj_value = float(line[len(obj_txt):])
            if args.solval:
                if abs(args.solval-obj_value) > 0.001:
                    ok=False
                    print("Wrong objective value, expected {} but got {}.".format(args.solval,obj_value))
            continue
        elif line.startswith('#'):
            continue
        v,val_ = line.strip().split()
        #print(v,val_)
        val=int(val_)
        if v not in valdict:
            continue
        if valdict[v]!=val:
            ok=False
            print("Wrong value, expexted {}={}, but got {}.".format(v,valdict[v],val))
        del valdict[v]
        
if len(valdict)>0:
    print("Warning, not all values were checked. Unchecked: {}".format(" ".join(valdict.keys())))
    ok=False

if not ok:
    sys.exit(1)
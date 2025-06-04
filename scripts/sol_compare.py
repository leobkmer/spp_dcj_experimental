from argparse import ArgumentParser
import sys
parser = ArgumentParser()

parser.add_argument("ding")
parser.add_argument("spp2")

args = parser.parse_args()

with open(args.ding) as df:
            for line in df:
                if not "Objective value" in line:
                    continue
                dingsol = int(line.strip().split("=")[1])
with open(args.spp2) as sf:
            for line in sf:
                if not "Objective value" in line:
                    continue
                sppsol = int(line.strip().split("=")[1])

print(dingsol,sppsol,file=sys.stderr)
if dingsol!=sppsol:
    print("Incongruent solution: ",dingsol,sppsol,file=sys.stderr)
    sys.exit(1)

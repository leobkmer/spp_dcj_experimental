from argparse import ArgumentParser


parser = ArgumentParser()

parser.add_argument("infile")

args = parser.parse_args()


with open(args.infile) as f:
    for line in f:
        if line.startswith(">"):
            line = '_'.join(line.split('_')[:-1])+"\n"
        print(line,end='')

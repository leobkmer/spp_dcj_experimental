#!/usr/bin/env python3

from sys import argv
import csv



lines = dict()
header = ''
for f in argv[1:]:
    for line in csv.reader(open(f), delimiter='\t'):
        if line[0].startswith('#'):
            header = '\t'.join(line)
        else:
            lines[tuple(line[:-1])] = line[-1]

print(header)

for line, val in lines.items():
    print('\t'.join(line + (val, )))

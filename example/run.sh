#!/usr/bin/env bash
python3 ../scripts/spp_dcj.py -a 1 -m idmap.txt tree.txt adjacencies.txt > example.ilp
python3 ../scripts/compute_initial.py -a 1 tree.txt adjacencies.txt idmap.txt example_init.sol
gurobi_cl InputFile=example_init.sol ResultFile=example.sol example.ilp 
python3 ../scripts/sol2adjacencies.py example.sol idmap.txt  > resolved_adjacencies.txt

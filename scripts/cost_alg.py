from ete3 import Tree
import pandas as pd
import argparse
from collections import Counter
import os

def initialise_gene_trees(tree, marker_set, counts):
    multi_count = lambda gene: True if sum([(counts[entry][gene]==1) for entry in counts.keys()]) != len(counts.keys()) else False
    accessory = [gene for gene in list(marker_set) if multi_count(gene)]
    all_trees = {gene:tree.copy() for gene in accessory}
    for gene in all_trees.keys(): 
        for leaf in all_trees[gene].iter_leaves():
            count = counts[leaf.name][gene]
            states = (count, count)
            leaf.add_feature("states", states)
    return [gene for gene in list(marker_set) if not gene in accessory], all_trees

def read_unimog(unimog_path, tree):
    genomes = {}
    marker_set = set()
    strip_orient = lambda x: x.replace("-", "").replace("+", "")
    with open(unimog_path, "r") as unimog:
        for line in unimog:
            if line.startswith('>'):
                entry = line[1:].strip()
            else:
                genome = line.strip(" |)\n").split(" ")   
                if entry in genomes.keys():   
                    genomes[entry].extend(list(map(strip_orient, genome)))
                else:
                    genomes[entry] = list(map(strip_orient, genome))
                marker_set.update(genomes[entry]) 
    counts = {entry:Counter(genomes[entry]) for entry in genomes.keys()}
    return initialise_gene_trees(tree, marker_set, counts)

def read_zombi(zombi_path, tree):
    marker_set = set()
    counts = {}
    for leaf in tree.iter_leaves():
        genome_df = pd.read_csv(f"{zombi_path}/Genomes/{leaf.name}_GENOME.tsv", sep="\t")
        counts[leaf.name] = genome_df["GENE_FAMILY"].value_counts()
        marker_set.update(counts.index)
    return initialise_gene_trees(tree, marker_set, counts)

def read_tree(filepath): #read in Newick tree with unlabelled internal nodes, and label the internal nodes
    t = Tree(filepath, format=1)
    i=1
    for n in t.traverse():
        if not n.is_leaf():
            n.name = f"I_{i}"
            i+=1
    return t

def recon_bottom_up(node, cost):
    if node.is_leaf():
        node.add_feature("possible_states", node.states)  # Assign leaf state
    else:
        left_states = recon_bottom_up(node.children[0], cost)
        right_states = recon_bottom_up(node.children[1], cost)

        # Ensure the possible_states attribute exists
        node.add_feature("possible_states", dict())

        if left_states[0]<=right_states[0]:
            node.possible_states = (min(min(right_states[1], left_states[1]), right_states[0]), max(min(right_states[1], left_states[1]), right_states[0]))
            cost[0] += max(0,right_states[0]-left_states[1])
        else:
            node.possible_states = (min(left_states[0], min(right_states[1], left_states[1])), max(left_states[0], min(right_states[1], left_states[1])))
            cost[0] += max(0,left_states[0]-right_states[1])

    return node.possible_states

def recon_top_down(node, parent=None):

    if node.is_leaf():
        node.add_feature("final_states", node.possible_states)
    elif parent is None:
        node.add_feature("final_states", node.possible_states)
    else:
        node.add_feature("final_states", dict())
        # Step I: Check if preliminary nodal set ⊆ parent final set, i.e. check if parent was formed through intersection
        if parent.final_states[0]>=node.possible_states[0] and parent.final_states[1]<=node.possible_states[1]:
            # Step II: Eliminate markers not in parent’s final set
            node.final_states = parent.final_states
        else:
            # Step III: Determine whether it was a union or intersection
            left_states = node.children[0].possible_states
            right_states = node.children[1].possible_states
            
            # Step IV: add parental states shared with at least one child
            if left_states[0]<=right_states[0]: #determine which child is to the left or right of node interval
                final_set = (min(max(left_states[0], parent.final_states[0]), node.possible_states[0]), max(min(right_states[1], parent.final_states[1]),node.possible_states[0]))
            else:
                final_set = (min(max(right_states[0], parent.final_states[0]), node.possible_states[0]), max(min(left_states[1], parent.final_states[1]),node.possible_states[1]))
            node.final_states = final_set

    # Recursively process child nodes
    for child in node.children:
        recon_top_down(child, node)

def reconstruction(tree):
    """Executes both phases on an ete3 Tree."""
    root = tree.get_tree_root()
    score = [0]
    recon_bottom_up(root, score)
    recon_top_down(root)
    return score

def intervals_to_tsv(t, recon_trees, core, tsv):
    with open(tsv, "w") as f:
        f.write("genome\tfamily\tlower\thigher\n")
        for block in recon_trees.keys():
            tree = recon_trees[block]
            for node in tree.traverse():
                if not node.is_leaf():
                    counts = (tree&node.name).final_states
                    lower = counts[0]
                    upper = counts[1]
                    if not (lower==0 and upper==0):
                        f.write(f"{node.name}\t{block}\t{lower}\t{upper}\n")
        for node in t.traverse():
            if not node.is_leaf():
                for block in core:
                    f.write(f"{node.name}\t{block}\t1\t1\n")

def output_scores(scores, tsv):
    with open(tsv, "w") as f:
        f.write("block\tscore\n")
        for block in scores.keys():
            f.write(f"{block}\t{scores[block][0]}\n")

def print_tree(node, level=0):
    """Utility function to print the tree with assigned states."""
    print("  " * level + f"{node.name}: {node.final_states}, {node.possible_states}")
    for child in node.children:
        print_tree(child, level + 1)

parser = argparse.ArgumentParser()
parser.add_argument('genomes', help="path to unimog file of genomes or zombi output directory")
parser.add_argument('tree', help="tree in newick format")
parser.add_argument('--input_format', choices=['zombi', 'unimog'], default='unimog')
parser.add_argument('--output_dir', help="directory to write output to", default=os.getcwd())

args = parser.parse_args()

tree = Tree(args.tree, format=1)

if args.input_format=="unimog":
    core, all_trees = read_unimog(args.genomes,tree)
else:
    core, all_trees = read_zombi(args.genomes,tree)

scores = {}
for block in all_trees.keys():
    scores[block] = reconstruction(all_trees[block])
    #print_tree(all_trees[block])
    #print(scores[block])

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)
intervals_to_tsv(tree, all_trees,core,f"{output_dir}/ranges.tsv")
output_scores(scores, f"{output_dir}/scores.tsv")
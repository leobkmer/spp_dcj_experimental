from argparse import ArgumentParser
import sys
def read_adjacencies(fl):
    adj = dict()
    with open(fl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomea, markera, extremitya, _, markerb, extremityb, _ =  line.strip().split()
            if not genomea in adj:
                adj[genomea]=[]
            adj[genomea].append(((markera,extremitya),(markerb,extremityb)))
        return adj

def get_singular_can_adjacencies(adj):
    singular = dict()
    for genome, adjacencies in adj.items():
        singular[genome]=set()
        marker_per_fam = dict()
        for (markera,extremitya),(markerb,extremityb) in adjacencies:
            fa = markera.split("_")[0]
            fb = markerb.split("_")[0]
            if not fa in marker_per_fam:
                marker_per_fam[fa]=set()
            if not fb in marker_per_fam:
                marker_per_fam[fb]=set()
            marker_per_fam[fa].add(markera)
            marker_per_fam[fb].add(markerb)
        for (markera,extremitya),(markerb,extremityb) in adjacencies:
            fa = markera.split("_")[0]
            fb = markerb.split("_")[0]
            canonical = tuple(sorted(((fa,extremitya),(fb,extremityb))))
            if len(marker_per_fam[fa])>1 or len(marker_per_fam[fb])>1:
                continue
            singular[genome].add(canonical)
    return singular


def get_leaves(pctable):
    all_nodes = set()
    not_leaves = set()
    with open(pctable) as f:
        for line in f:
            child,parent = line.strip().split()
            all_nodes.add(child)
            all_nodes.add(parent)
            not_leaves.add(parent)
        leaves = all_nodes.difference(not_leaves)
    return leaves
            
        

    


parser = ArgumentParser()
parser.add_argument("truth")
parser.add_argument("prediction")
parser.add_argument("tree")

args = parser.parse_args()

leaves = get_leaves(args.tree)

truth_a = get_singular_can_adjacencies(read_adjacencies(args.truth))
pred_a = get_singular_can_adjacencies(read_adjacencies(args.prediction))

for genome, padjset in pred_a.items():
    if genome in leaves:
        continue
    if not genome in truth_a:
        print("Warning: genome",genome, "not found in ground truth. Skipping",file=sys.stderr)
        continue
    tadjset = truth_a[genome]
    true_positives = len(tadjset.intersection(padjset))
    false_positives = len(padjset.difference(tadjset))
    false_negatives = len(tadjset.difference(padjset))

    print(genome,true_positives,false_positives,false_negatives)
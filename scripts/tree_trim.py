import data_utils as du




def invert_mrk(mrk):
    r = []
    for s,m in mrk:
        r.append(("+" if s=="-" else "-",m))
    return r[::-1]

def canonicalize_chrm(chrm):
    lc,mrk = chrm
    if lc == du.CHR_LINEAR:
        return (du.CHR_LINEAR,min(mrk,invert_mrk(mrk)))
    assert(lc==du.CHR_CIRCULAR)
    #circular chromosome
    #inefficient, but for simplicity's sake
    rotations=[mrk[i:]+mrk[:i] for i in range(len(mrk))]
    rrotatiions = [invert_mrk(y) for y in rotations]
    return (du.CHR_CIRCULAR,min(rotations+rrotatiions))

def canonicalize_gnm(gnm):
    name,chrs = gnm
    can_chrs = [canonicalize_chrm(chrm) for chrm in chrs]
    return (name,sorted(can_chrs,key=lambda c: (sorted([m for _,m in c[1]]),c)))





def canon_adjacencies(gnm):
    return [tuple(sorted(adj)) for adj in du.unimog2adjacencies(canonicalize_gnm(gnm))]


def canonicalized_has_distance0or1(adja,adjb):
    #check if we find a distance of 0 or 1 based on canonical adjacency representations
    #Note: not all distance 1 cases can be found this way
    diffa = set(adja).difference(set(adjb))
    diffb = set(adjb).difference(set(adja))
    #print(diffa)
    #print(diffb)
    if len(diffa)==0 and len(diffb)==0:
        return 0
    if len(diffa)==1 and len(diffb)==2:
        axs=list(sorted(diffa.pop()))
        bxs = list(sorted([(m,x) for y in diffb for (m,x) in y if x!= du.EXTR_CAP]))
        if axs==bxs:
            return 1
    if len(diffb)==1 and len(diffa)==2:
        bxs=list(sorted(diffb.pop()))
        axs = list(sorted([(m,x) for y in diffa for (m,x) in y if x!= du.EXTR_CAP]))
        if axs==bxs:
            return 1
        

    if len(diffa)!=2:
        return None
    if len(diffb)!=2:
        return None
    ((a,b),(c,d))=tuple(sorted([sorted(x) for x in diffa]))
    ((e,f),(g,h))=tuple(sorted([sorted(x) for x in diffb]))
    if not a==e:
        return None
    if c==f and tuple(sorted((b,d)))==(g,h):
        return 1
    if d==f and tuple(sorted((b,c)))==(g,h):
        return 1
    return None

UNIVERSE = "all"
def trim_tree_fitch(tree,leaves,unimog):
    filtered_nodes = []
    root,pc_tree=du.cp_to_pc(tree)
    possible_linearizations = dict()
    for gnm in unimog:
        name,_=gnm
        ca = canon_adjacencies(gnm)
        possible_linearizations[name]=[ca]
    traversal = bottom_up_traversal(root,pc_tree)
    for name in traversal:
        if name in leaves:
            continue
        children = pc_tree[name]
        if len(children) > 2:
            possible_linearizations[name]=UNIVERSE
        if len(children)==1:
            possible_linearizations[name]=possible_linearizations[children[0]]
        if len(children)==2:
            c1 = children[0]
            c2 = children[1]
            if possible_linearizations[c1]==UNIVERSE or possible_linearizations[c2]==UNIVERSE:
                possible_linearizations[name]=UNIVERSE
            else:
                fitchable = True
                has_intersection = False
                for a in possible_linearizations[c1]:
                    for b in possible_linearizations[c2]:
                        d = canonicalized_has_distance0or1(a,b)
                        if d is None:
                            fitchable=False
                        if d==0:
                            has_intersection=True
                if not fitchable:
                    possible_linearizations[name]=UNIVERSE
                elif has_intersection:
                    filtered_nodes.append(c1)
                    filtered_nodes.append(c2)
                    possible_linearizations[name]=[x for x in possible_linearizations[c1] if x in possible_linearizations[c2]]
                else:
                    filtered_nodes.append(c1)
                    filtered_nodes.append(c2)
                    possible_linearizations[name]=possible_linearizations[c1]+possible_linearizations[c2]
    return filtered_nodes,possible_linearizations

def bottom_up_traversal(root, pc_tree):
    traversal = [root]
    curr_level = [root]
    while True:
        next_level = []
        for x in curr_level:
            next_level.extend(pc_tree.get(x,[]))
        traversal.extend(next_level)
        curr_level=next_level
        if len(curr_level)==0:
            break
    #bottom up traversal
    traversal=traversal[::-1]
    return traversal
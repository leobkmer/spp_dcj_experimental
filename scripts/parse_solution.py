from argparse import ArgumentParser
import data_utils as du
import sys
import networkx as nx


EPSILON = 0.0001
def read_solfile(solfile):
    var_val = {}
    with open(solfile) as sf:
        for line in sf:
            if line.startswith("#"):
                #comment
                continue
            varname,val = line.strip().split()
            var_val[varname]=float(val)
    return var_val
                
                
def annotate_diagram(relationalDiagram,var_val,tree_ids):
    for te,G in relationalDiagram["graphs"].items():
        tei = tree_ids[te]
        for u,v,k,data in G.edges(keys=True,data=True):
            xvar = "x{sep}{tei}{sep}{eid}".format(sep=du.SEP,tei=tei,eid=data['id'])
            if xvar not in var_val:
                print(G.nodes[u]['type'],G.nodes[v]['type'])
                assert(G.nodes[u]['type']==du.VTYPE_CAP or G.nodes[v]['type']==du.VTYPE_CAP)
                assert(data['type']==du.ETYPE_ADJ)
                print("Warning: Telomere edge {} does not appear in solution. Defaulting to 0.".format(xvar),file=sys.stderr)
                G[u][v][k]["x"]=0
                if data['type']==du.ETYPE_ADJ:
                    G[u][v][k]["a"]=0

            else:
                G[u][v][k]["x"]=var_val[xvar]
                if data['type']==du.ETYPE_ADJ:
                    avar = "a{sep}{anc}".format(sep=du.SEP,tei=tei,anc=data['id'])
                    G[u][v][k]["a"]=var_val[avar]

        for u,data in G.nodes(data=True):
            gvar = "g{sep}{anc}".format(sep=du.SEP,anc=data['anc'])
            if gvar not in var_val:
                assert(G.nodes[u]['type']==du.VTYPE_CAP)
                print("Warning: Telomere {u} does not appear in solution. Defaulting to 0.".format(u=u),file=sys.stderr)
                G.nodes[u]["g"]=0
            else:
                G.nodes[u]["g"]=var_val[gvar]


                
def sanity_check_decomposition(relationalDiagram,var_val,tree_ids,fam_bounds,famsep='_',check_fam_bounds=True):
    annotate_diagram(relationalDiagram,var_val,tree_ids)
    for te,G in relationalDiagram["graphs"].items():
        for u,v,k,data in G.edges(keys=True,data=True):
            assert(data["x"] - G.nodes[u]["g"] <= EPSILON)
            assert(data["x"] - G.nodes[v]["g"] <= EPSILON)
            assert(data["x"] >= -EPSILON and data["x"] <= 1+EPSILON)
            if data['type']==du.ETYPE_ADJ:
                assert(data["a"] >= -EPSILON and data["a"] <= 1+EPSILON)
                assert(abs(data["x"]-data["a"])<=EPSILON)
        ga,gb = te
        xtr_multiplicities = {ga : {},gb:{}}
        for u in G.nodes():
            assert(G.nodes[u]["g"] >=-EPSILON and G.nodes[u]["g"] <= 1+EPSILON)
            acands = [v for v in G[u] for k in G[u][v] if abs(1-G[u][v][k]["x"])<=EPSILON and G[u][v][k]["type"]==du.ETYPE_ADJ]
            ecands = [v for v in G[u] for k in G[u][v] if abs(1-G[u][v][k]["x"])<=EPSILON and G[u][v][k]["type"] in [du.ETYPE_EXTR,du.ETYPE_ID]]
            if abs(1-G.nodes[u]["g"])<=EPSILON:
                assert(len(acands)==1)
                if G.nodes[u]['type']!=du.VTYPE_CAP:
                    assert(len(ecands)==1)
            else:
                assert(len(acands)==0)
                assert(len(ecands)==0)
            if G.nodes[u]['type']==du.VTYPE_EXTR and abs(1-G.nodes[u]["g"])<=EPSILON:
                (gnm,(gene,_))=G.nodes[u]['id']
                fam=gene.split(famsep)[0]
                if not fam in xtr_multiplicities[gnm]:
                    xtr_multiplicities[gnm][fam]=0
                xtr_multiplicities[gnm][fam]+=1
        if check_fam_bounds:
            for gnm,mlts in xtr_multiplicities.items():
                for fm,n in mlts.items():
                    lo,hi = fam_bounds[gnm][fm]
                    assert(n%2==0)
                    assert(lo <= int(n/2))
                    assert(hi >= int(n/2))
        for eid,did in relationalDiagram['siblings'][te]:
            tei = tree_ids[te]
            ex = "x{sep}{te}{sep}{eid}".format(eid=eid,sep=du.SEP,te=tei)
            dx = "x{sep}{te}{sep}{did}".format(did=did,sep=du.SEP,te=tei)
            assert(abs(var_val[ex]-var_val[dx])<=EPSILON)
    

    
        


def make_edge_id_to_edge_maps(relationalDiagrams):
    emaps = dict()
    for e,G in relationalDiagrams['graphs']:
        emaps[e]=dict()
        for u,v,k,data in G.edges(keys=True,data=True):
            emaps[e][data['id']]=(u,v,k)
    return emaps

def add_standard_arguments(parser):
    parser.add_argument('phylogeny_edge_ids',
            help='phylogenetic tree as child->parent relation table with their assigned ids by spp_dcj.py (Important: children must be in first column.)')
    parser.add_argument('candidateAdjacencies',
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('idmap',
            help='The table with ID-to-gene extremity mapping')
    parser.add_argument("solfile")
    parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
    grp= parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('-t', '--no_telomeres', action='store_true',
            help='don\'t add any additional telomeres. Overrides -at.')
    grp.add_argument('-l', '--all_telomeres', action='store_true',
            help='Add additional telomeres to all nodes.')
    grp.add_argument('-d',action='store_true',help="Default mode.")
    grp.add_argument('--guess',action='store_true',help="*Only use if you know what you're doing.* Try to guess added telomeres from id file.")

def guess_telomere_adding(glob_table):
    gnmwise = dict()
    for (gnm,(gne,xtr)) in glob_table.keys():
        if not gnm in gnmwise:
            gnmwise[gnm]={"cap":set(),"nocap":set()}
        if xtr == du.EXTR_CAP:
            gnmwise[gnm]["cap"].add(gne)
        else:
            gnmwise[gnm]["nocap"].add(gne)
    t_dont_add=True
    t_add_all=False
    for gnm,d in gnmwise.items():
        n_added = 0
        for x in d["cap"]:
            assert(x.startswith("t_"))
            suffix = x[2::]
            if suffix in d["nocap"]:
                n_added+=1
        if n_added >= 1:
            t_dont_add=False
        elif n_added == len(d["nocap"]):
            t_add_all=True
    if t_dont_add:
        print("Guessing has identified file as having no added telomeres",file=sys.stderr)
    if t_add_all:
        print("Guessing has identified file as having all added telomeres",file=sys.stderr)
    if not t_dont_add and not t_add_all:
        print("Guessing has identified file as having default mode",file=sys.stderr)
    return t_dont_add,t_add_all


def sanity_check_groups(homology_graph,fam_sep):
    for c in nx.connected_components(homology_graph):
        observed_families = set()
        base_name = None
        for gnm,gne in c:
            assert(gnm not in observed_families)
            observed_families.add(gnm)
            curr_basename = gne.split(fam_sep)[0]
            if base_name is None:
                base_name=curr_basename
            assert(base_name == curr_basename)

            

def get_relabeled_group_ids(relationalDiagrams,fam_sep='_'):
    new_group_ids = dict()
    homology_graph = nx.Graph()
    for te,G in relationalDiagrams["graphs"].items():
        for _,data in G.nodes(data=True):
            if data["type"]!=du.VTYPE_CAP and abs(1-data["g"])>=0:
                (ugnm,(ugne,_)) = data["id"]
                homology_graph.add_node((ugnm,ugne))
        for u,v,k,data in G.edges(keys=True,data=True):
            if data['type']==du.ETYPE_EXTR and abs(1-data['x']) <= EPSILON:
                if G.nodes[u]["type"]!=du.VTYPE_CAP:
                    (ugnm,(ugne,_)) = G.nodes[u]["id"]
                    (vgnm,(vgne,_)) = G.nodes[v]["id"]
                    homology_graph.add_edge((ugnm,ugne),(vgnm,vgne))
    sanity_check_groups(homology_graph,fam_sep)
    fam_counts = dict()
    for c in nx.connected_components(homology_graph):
        c = list(c)
        basename= c[0][1].split(fam_sep)[0]
        if not basename in fam_counts:
            fam_counts[basename]=0
        fam_counts[basename]+=1
        newname = "{basename}{famsep}{mx}{famcount}".format(basename=basename,famsep=fam_sep,famcount=fam_counts[basename],mx="m" if len(c)>=2 else "x")
        for x in c:
            new_group_ids[x]=newname
    return new_group_ids
        
        

def get_adjacencies(relationalDiagrams):
    adj = set()
    for te,G in relationalDiagrams["graphs"].items():
        for u,v,data in G.edges(data=True):
            if data["type"]!=du.ETYPE_ADJ:
                continue
            if abs(1-data["a"]) <= EPSILON:
                u_=G.nodes[u]["id"]
                v_=G.nodes[v]["id"]
                u_,v_ = (u_,v_) if u_ < v_ else (v_,u_)
                adj.add((u_,v_))
    return adj
            
    
def rename_adjacencies(adj,new_ids):
    new_adj = []
    for (agnm,(agne,axtr)),(bgnm,(bgne,bxtr)) in adj:
        nagne = new_ids[(agnm,agne)] if axtr != du.EXTR_CAP else agne
        nbgne = new_ids[(bgnm,bgne)] if bxtr != du.EXTR_CAP else bgne
        new_adj.append(((agnm,(nagne,axtr)),(bgnm,(nbgne,bxtr))))
    return new_adj
        
    




        


    

def main():
    parser = ArgumentParser(description="Parse a solution and write result files. If no result file is given, basic sanity checks are performed.")
    add_standard_arguments(parser)
    parser.add_argument('--family-bounds','-fmb')
    parser.add_argument('--write-adjacencies')
    parser.add_argument('--write-unimog')
    parser.add_argument('--write-resolved-families')
    parser.add_argument('--no-relabel-adjacencies',action='store_true')
    args = parser.parse_args()
    
    tree_ids_r = du.read_tree_edge_name_map(args.phylogeny_edge_ids)
    tree_ids = dict()
    for k,v in tree_ids_r.items():
        tree_ids[v]=k
    with open(args.candidateAdjacencies) as caf:
        candidateAdjacencies = du.parseAdjacencies(caf,sep=args.separator)
    with open(args.idmap) as idm:
        glob_table, loc_tables = du.parse_id_file(idm)
        global_ext2id = du.IdManager(0,is_cap=lambda x: False,initial_table=glob_table)
    # add telomeres
    tree_edges = list(tree_ids_r.values())
    leaves = set([x for x, v in du.getLeaves(tree_edges).items() if v])
    t_dont_add = False
    t_add_all = False
    if args.no_telomeres:
        t_dont_add=True
    elif args.all_telomeres:
        t_add_all=True
    elif args.guess:
        t_dont_add,t_add_all = guess_telomere_adding(glob_table)

    telomeres = du.identifyCandidateTelomeres(candidateAdjacencies,0,leaves,dont_add=t_dont_add,addToAll=t_add_all)
    # construct adjacency graphs
    genes = candidateAdjacencies['genes']
    adjacencies = candidateAdjacencies['adjacencies']
    weights = candidateAdjacencies['weights']
    fam_bounds = dict()
    if args.family_bounds:
        fam_bounds=du.parseFamilyBounds(args.family_bounds)
    families = candidateAdjacencies['families']
    du.fillFamilyBounds(families,fam_bounds)
    relationalDiagrams = du.constructRelationalDiagrams(tree_edges,
            adjacencies, telomeres, weights, genes, global_ext2id,fam_bounds=fam_bounds,
            sep=args.separator,loc_manager_tables=loc_tables)
    for G in relationalDiagrams['graphs'].values():
        du.checkGraph(G,cf=True)
    var_val = read_solfile(args.solfile)
    print("Checking solfile...",file=sys.stderr)
    sanity_check_decomposition(relationalDiagrams,var_val,tree_ids,fam_bounds,famsep=args.separator,check_fam_bounds=True if args.family_bounds else False)    
    new_ids = get_relabeled_group_ids(relationalDiagrams,fam_sep=args.separator)
    adjs_raw = get_adjacencies(relationalDiagrams)
    new_adjs = rename_adjacencies(adjs_raw,new_ids)
    print("Solfile is valid.",file=sys.stderr)
    if args.write_adjacencies:
        if args.no_relabel_adjacencies:
            write_adjs = adjs_raw
        else:
            write_adjs = new_adjs
        species_adj, weights = du.transform_list_toadj_dicts(write_adjs)
        with open(args.write_adjacencies,"w") as fl:
            du.writeAdjacencies(species_adj,weights,fl)
    if args.write_resolved_families:
        with open(args.write_resolved_families,"w") as fl:
            du.write_res_families(new_ids, fl)
    if args.write_unimog:
        genomes = du.canonicize_unimog(du.matched_adjacencies_to_unimog(new_adjs))
        with open(args.write_unimog,"w") as fl:
            du.write_unimog(genomes,fl)
    
    


    
    




if __name__ == "__main__":
    main()
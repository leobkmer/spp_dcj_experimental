import data_utils as du
import csv

def annotate_external_warm_start(leaves,graphs,adjacencies,matchings):
    print(adjacencies,matchings)
    for (a,b),g in graphs.items():
        for u,v,k,etype in g.edges(keys=True,data='type'):
            g[u][v][k]['is_set']=False
        #TODO: set edges based on adjacencies and matching
            (gnmu,(gu,xu)) = g.nodes[u]['id']
            (gnmv,(gv,xv)) = g.nodes[v]['id']

            if etype==du.ETYPE_ADJ:
                g[u][v][k]['is_set'] = ((gu,xu),(gv,xv)) in adjacencies.get(gnmu,set()) or ((gv,xv),(gu,xu)) in adjacencies.get(gnmu,set()) or gnmu in leaves
            
            if etype==du.ETYPE_EXTR:
                if gnmv!=a:
                    gu_,gv_= gv,gu
                else:
                    gu_,gv_ = gu,gv
                g[u][v][k]['is_set']=  (gu_,gv_) in matchings[(a,b)]
        #set indel edges of markers that were not matched
        for u,v,k,etype in g.edges(keys=True,data='type'):
            if etype != du.ETYPE_ID:
                continue
            mcands = [w for w in g[u] for l in g[u][w] if g[u][w][l]['etype']==du.ETYPE_EXTR and g[u][w][l]['is_set']]
            g[u][v][k]['is_set']=len(mcands)==0


#tab-separated format: genome,marker,xtr,marker,xtr
def parse_adjacencies(af):
    adjacencies = dict()
    with open(af) as f:
        reader = csv.reader(f,delimiter='\t')
        for entries in reader:
            print(entries)
            gnm,m1,x1,m2,x2 = entries
            if not gnm in adjacencies:
                adjacencies[gnm] = set()
            adjacencies[gnm].add(((m1,x1),(m2,x2)))
    return adjacencies

#tab-separated format: genome,marker,genome,marker
def parse_matchings(mf):
    matchings = dict()
    with open(mf) as f:
        for entries in csv.reader(f,delimiter='\t'):
            gnm1,m1,gnm2,m2 = entries
            if (gnm1,gnm2) not in matchings:
                matchings[(gnm1,gnm2)]=set()
                matchings[(gnm2,gnm1)]=set()
            matchings[(gnm1,gnm2)].add((m1,m2))
            matchings[(gnm2,gnm1)].add((m2,m1))
    return matchings
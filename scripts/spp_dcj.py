#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from sys import stdout, stderr, exit
from itertools import product, combinations, chain, repeat
from functools import reduce
from collections import defaultdict
from math import ceil
import logging
import csv

# import from third-party packages
import networkx as nx

# import from own packages
import data_utils as du


#
# global variables
#

ADJ_TRUST_THRESHOLD = 0.9

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

#
# ILP OBJECTIVE
#

def objective(graphs, circ_singletons, alpha, out):
    out.write('maximize ')
    terms = ["{ialph} w{sep}{te} + {alph} f{sep}{te}".format(
        alph=alpha,ialph=1-alpha,sep=du.SEP,te=tree_edge)
        for tree_edge, _ in enumerate(sorted(graphs.items()))]
    print(" + ".join(terms),file=out)
        


#
# ILP CONSTRAINTS
#

def constraints(graphs, siblings, out):
    out.write('subject to\n')
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing constraints for relational diagram of {} and ' + \
                '{}').format(child, parent))
        genomes = [child,parent]
        global_constraints(G,out)
        loose_constraints(G,i,out)
        slm_constraints(G,i,siblings,out)
        regular_reporting(G,i,genomes,out)
        pcap_reporting(G,i,genomes,out)
        cs_constraints(G,i,out,max_circ_len=len(G.nodes())/2)

    out.write('\n')

def get_gene_extremities(G):
    ids_to_v = {data['id'] : v for v,data in G.nodes(data=True)}
    l = []
    for v,data in G.nodes(data=True):
        if data['type'] == du.VTYPE_CAP:
            continue
        l.append((v,ids_to_v[du.complement_id(data['id'])]))
    return l


#TODO: Integrate multiplicities!
def c01(G,out):
    for v,data in G.nodes(data=True):
        print("g{sep}{v} = 1".format(sep=du.SEP,v=v),file=out)


def c02(G,out):
    for u,v in get_gene_extremities(G):
        print("g{sep}{u} - g{sep}{v} = 0".format(u=u,v=v,sep=du.SEP),file=out)

def c03(G,out):
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            continue
        aes = [G[v][u][i]['id'] for u in G[v] for i in G[v][u] if G[v][u][i]['type']==du.ETYPE_ADJ]
        sm = ' + '.join(['a{sep}{e}'.format(sep=du.SEP,e=e) for e in aes])
        print("{sm} - g{sep}{v}".format(sm=sm,sep=du.SEP,v=v),file=out)


def global_constraints(G,out):
    for c in [c01,c02,c03]:
        c(G,out)

def c04(G,tree_edge,out):
    ws = [ '{w} x{sep}{te}{sep}{e}'.format(q=data['weight'],te=tree_edge,e=data['id'],sep=du.SEP) for u,v,data in G.edges(data=True) if data['type']==du.ETYPE_ADJ]
    print("{s} - w{sep}{treeedge} = 0".format(treeedge=tree_edge,s=' + '.join(ws),sep=du.SEP),file=out)


def c05(G,tree_edge,out):
    print("n{sep}{te} - c{sep}{te} + q{sep}{te} + s{sep}{te} - f{sep}{te} = 0".format(sep=du.SEP,te=tree_edge),file=out)

def c06(G,tree_edge,out):
    xs = ['0.5 x{sep}{te}{sep}{e}'.format(te=tree_edge,sep=du.SEP,e=data['id']) for u,v,data in G.edges(data=True) if data['type']==du.ETYPE_ADJ]
    print('{s} - n{sep}{te} = 0'.format(s=' + '.join(xs),sep=du.SEP,te=tree_edge),file=out)

def c07(G,tree_edge,out):
    cs = ['rc{sep}{te}{sep}{v}'.format(te=tree_edge,sep=du.SEP,v=v) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP]
    print('{s} - c{sep}{te} = 0'.format(s=' + '.join(cs),sep=du.SEP,te=tree_edge),file=out)

def c08(G,tree_edge,out):
    print("pab{sep}{te} + pABa{sep}{te} + pABb{sep}{te} - pAB{sep}{te} - 2q{sep}{te}".format(sep=du.SEP,te=tree_edge),file=out)


def summary_constraint(pathtype,predicate,G,tree_edge,out):
    sm = ['r{ptype}{sep}{te}{sep}{v}'.format(te=tree_edge,ptype=pathtype,sep=du.SEP,v=data['type']) for v,data in G.nodes(data=True) if predicate(data['type'])]
    print("{s} - p{ptype}{sep}{te}".format(s=' + '.join(sm),ptype=pathtype,sep=du.SEP,te=tree_edge),file=out)
    
c09 = lambda G,tree_edge,out: summary_constraint('ab',lambda x: x!=du.VTYPE_CAP,G,tree_edge,out)
c10 = lambda G,tree_edge,out: summary_constraint('Ab',lambda x: x==du.VTYPE_CAP,G,tree_edge,out)    
c11 = lambda G,tree_edge,out: summary_constraint('Ba',lambda x: x==du.VTYPE_CAP,G,tree_edge,out)    
c12 = lambda G,tree_edge,out: summary_constraint('Aa',lambda x: x==du.VTYPE_CAP,G,tree_edge,out)    
c13 = lambda G,tree_edge,out: summary_constraint('Bb',lambda x: x==du.VTYPE_CAP,G,tree_edge,out)    

def c14toc17(G,tree_edge,out):
    for indel in "ab":
        for tel in "AB":
            print("pAB{ind}{sep}{te} - p{tel}{ind}{sep}{te}".format(ind=indel,tel=tel,te=tree_edge,sep=du.SEP),file=out)
c18 = lambda G,tree_edge,out: summary_constraint('AB',lambda x: x==du.VTYPE_CAP,G,tree_edge,out)

def c19(G,tree_edge,out):
    sm = ['rs{sep}{te}{sep}{v}'.format(te=tree_edge,sep=du.SEP,v=data['type']) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP]
    print("{s} - s{sep}{te}".format(s=' + '.join(sm),sep=du.SEP,te=tree_edge),file=out)

def c20(G,tree_edge,out):
    for v in G.nodes():
        ees = [G[v][u][i]['id'] for u in G[v] for i in G[v][u] if G[v][u][i]['type'] in [du.ETYPE_EXTR,du.ETYPE_ID]]
        sm = ' + '.join(['x{sep}{te}{sep}{e}'.format(sep=du.SEP,te=tree_edge,e=e) for e in ees])
        print("{sm} - g{sep}{v} = 0".format(sm=sm,sep=du.SEP,v=v),file=out)

def c21(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if data['type'] == du.ETYPE_ADJ:
            print("a{sep}{e} - x{sep}{te}{sep}{e} = 0".format(sep=du.SEP,te=tree_edge,e=data['id']),file=out)

def c22(G,tree_edge,out):
    for v in G.nodes():
        print("z{sep}{te}{sep}{v} - g{sep}{v} <= 0".format(sep=du.SEP,te=tree_edge,v=v),file=out)


def loose_constraints(G,tree_edge,out):
    for c in [c04,c05,c06,c07,c08,c09,c10,c11,c12,c13,c14toc17,c18,c19,c20,c21,c22]:
        c(G,tree_edge,out)

def c23(sibs,tree_edge,out):
    for eid,did in sibs:
        print("x{sep}{te}{sep}{eid} - x{sep}{te}{sep}{did} = 0".format(sep=du.SEP,eid=eid,did=did,te=tree_edge),file=out)

def c24(G,tree_edge,out):
    for e,data in G.edges(data=True):
        eid = data['id']
        for u,v in [e,e[::-1]]:
            print("y{sep}{te}{sep}{v} - {u} x{sep}{te}{sep}{eid} - y{sep}{te}{sep}{u} >= {u}".format(sep=du.SEP,eid=eid,v=v,u=u,te=tree_edge),file=out)

def c25(G,tree_edge,out):
    for v in G.nodes():
        print("{v} z{sep}{te}{sep}{v} - y{sep}{te}{sep}{v} <= 0".format(sep=du.SEP,v=v,te=tree_edge),file=out)

def slm_constraints(G,tree_edge,sibs,out):
    """
    Shao-Lin-Moret Constraints as defined in Table 1.
    """
    c23(sibs,tree_edge,out)
    c24(G,tree_edge,out)
    c25(G,tree_edge,out)






def c26(G,tree_edge,genomes,out):
    for x,y, data in G.edges(data=True):
        if data['type'] == du.ETYPE_ID:
            for v in [x,y]:
                if get_genome(G, v) == genomes[0]:
                    print("l{sep}{te}{sep}{v} +  x{sep}{te}{sep}{e} <= 1".format(sep=du.SEP,te=tree_edge,v=v,e=data['id']),file=out)
                else:
                    print("l{sep}{te}{sep}{v} -  x{sep}{te}{sep}{e} >= 0".format(sep=du.SEP,te=tree_edge,v=v,e=data['id']),file=out)

def get_genome(G, v):
    return G.nodes[v]['id'][0]


def c27(G,tree_edge,genomes,out):
    for x,y, data in G.edges(data=True):
        if data['type'] == du.ETYPE_ID or du.VTYPE_CAP in [G.nodes[v]['type'] for v in [x,y]]:
            continue

        for u,v in [(x,y),(y,x)]:
            if get_genome(G,v) == genomes[0]  and data['type'] == du.ETYPE_ADJ:
                print("l{sep}{te}{sep}{v} - l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} - rab{sep}{te}{sep}{v} <= 1".format(
                    sep=du.SEP,
                    te=tree_edge,
                    e=data['id'],
                    u=u,
                    v=v
                ), file=out)
            else:
                print("l{sep}{te}{sep}{v} - l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                    sep=du.SEP,
                    te=tree_edge,
                    e=data['id'],
                    u=u,
                    v=v
                ), file=out)
            
def c28(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP or get_genome(G,v)!=genomes[0]:
            continue
        print('rc{sep}{te}{sep}{v} - z{sep}{te}{sep}{v} <= 0'.format(
            sep=du.SEP,
            te=tree_edge,
            v=v
        ), file=out)

def c29(G,tree_edge,genomes,out):
    for x, y, data in G.edges(data=True):
        if data[type] != du.ETYPE_ID or get_genome(G,x) != genomes[0]:
            continue
        for v in [x,y]:
            print("rab{sep}{te}{sep}{v} - x{sep}{te}{sep}{e} <= 0".format(
                sep=du.SEP,
                te=tree_edge,
                v=v,
                e=data['id']
            ),file=out)

def regular_reporting(G,tree_edge,genomes,out):
    for c in [c26,c27,c28,c29]:
        c(G,tree_edge,genomes,out)


def c30(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        print("l{sep}{te}{sep}{v} = {gid}".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            gid=0 if get_genome(G,v)==genomes[0] else 1
        ),file=out)

def c31(G,tree_edge,genomes,out):
    for x,y,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ADJ or not du.VTYPE_CAP in [G.nodes[v]['type'] for v in [x,y]]:
            continue
        if get_genome(G,x) == genomes[0]:
            if G.nodes[x]['type'] == du.VTYPE_CAP:
                v,u = x,y
            else:
                v,u = y,x
            print("l{sep}{te}{sep}{u} -  l{sep}{te}{sep}{v} - rAB{sep}{te}{sep}{v} - rAb{sep}{te}{sep}{v} + x{sep}{te}{sep}{e} <= 1".format(
                sep=du.SEP,
                te=tree_edge,
                u=u,
                v=v,
                e=data['id']
            ),file=out)
        else:
            if G.nodes[x]['type'] == du.VTYPE_CAP:
                u,v = x,y
            else:
                u,v = y,x
            print("l{sep}{te}{sep}{u} -  l{sep}{te}{sep}{v} - rBa{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                sep=du.SEP,
                te=tree_edge,
                u=u,
                v=v,
                e=data['id']
            ),file=out)

def c32(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if get_genome(G,v)==genomes[0] and data['type']==du.VTYPE_CAP:
            print("rAB{sep}{te}{sep}{v} - z{sep}{te}{sep}{v} <= 0".format(
                sep=du.SEP,
                v=v,
                te=tree_edge
            ),file=out)

def c33(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!=du.VTYPE_CAP:
            continue
        if get_genome(G,v)==genomes[0]:
            reps = ["Ab","Aa"]
        else:
            reps = ["Ba","Bb"]
        print("r{a}{sep}{te}{sep}{v} + r{b}{sep}{te}{sep}{v} + y{sep}{te}{sep}{v} >= 1".format(
            a=reps[0],
            b=reps[1],
            te=tree_edge,
            sep=du.SEP,
            v=v
        ),file=out)

def c34(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!=du.VTYPE_CAP:
            continue
        if get_genome(G,v)==genomes[0]:
            reps = ["Ab","Aa"]
        else:
            reps = ["Ba","Bb"]
        for r in reps:
            print("y{sep}{te}{sep}{v} + {v} r{r}{sep}{te}{sep}{v} <= {v}".format(
                sep=du.SEP,
                te=tree_edge,
                v=v,
                r=r
            ),file=out)

def c35(G,tree_edge,genomes,out):
    for x,y,data in G.edges(data=True):
        if data['type'] != du.ETYPE_ADJ:
            continue
        if G.nodes[x]['type'] == du.VTYPE_CAP:
            u,v = y,x
        elif G.nodes[y]['type'] == du.VTYPE_CAP:
            u,v = x,y
        else:
            continue
        if get_genome(G,v) == genomes[0]:
            for r in ["AB","Ab"]:
                print("r{r}{sep}{te}{sep}{v} -  l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                    r=r,
                    sep=du.SEP,
                    te=tree_edge,
                    v=v,
                    u=u,
                    e=data['id']
                ),file=out)
        else:
            print("rBa{sep}{te}{sep}{v} +  l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 2".format(
                    r=r,
                    sep=du.SEP,
                    te=tree_edge,
                    v=v,
                    u=u,
                    e=data['id']
                ),file=out)
            
def pcap_reporting(G,tree_edge,genomes,out):
    for c in [c30,c31,c32,c33,c34,c35]:
        c(G,tree_edge,genomes,out)


def c36(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if du.VTYPE_CAP in [G.nodes[x]['type'] for x in [u,v]]:
            continue
        if data['type']==du.ETYPE_EXTR:
            continue
        print("d{sep}{te}{sep}{v} + d{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 2".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            u=u,
            e=data['id']
        ),file=out)
        print("d{sep}{te}{sep}{v} + d{sep}{te}{sep}{u} - x{sep}{te}{sep}{e} >= 0".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            u=u,
            e=data['id']
        ),file=out)

def c37(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ID:
            continue
        print("w{sep}{te}{sep}{u} - w{sep}{te}{sep}{v} = 0".format(
            sep=du.SEP,
            te=tree_edge,
            u=u,
            v=v
        ),file=out)

def c38(G,tree_edge,out,max_circ_len):
    for u,v,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ADJ or du.VTYPE_CAP in [G.nodes[x]['type'] for x in [u,v]]:
            continue
        print("w{sep}{te}{sep}{u} - w{sep}{te}{sep}{v} + d{sep}{te}{sep}{v} - d{sep}{te}{sep}{u} + {K} x{sep}{te}{sep}{e} - {K} rs{sep}{te}{sep}{u} - {K} rs{sep}{te}{sep}{v} <= {K}".format(
            sep=du.SEP,
            u=u,
            v=v,
            e=data['id'],
            K=max_circ_len,
            te=tree_edge
        ),file=out)

def cs_constraints(G,tree_edge,out,max_circ_len):
    c36(G,tree_edge,out)
    c37(G,tree_edge,out)
    c38(G,tree_edge,out,max_circ_len)

'''
def c01(G, out):
    for v, vdata in G.nodes(data = True):
        line = ''
        for u in G.neighbors(v):
            # two vertices may share multiple edges
            for data in G[u][v].values():
                if data['type'] == du.ETYPE_ADJ:
                    line += line and ' + '
                    line += 'x{}'.format(data['id'])
        if line:
            if vdata['type'] == du.VTYPE_EXTR:
                line += ' = 1\n'
            elif vdata['type'] == du.VTYPE_CAP:
                line += ' - o{} = 0\n'.format(v)
            else:
                raise Exception('unknown node type!')

            out.write(line)


def c02(G, i, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for v in G.nodes():
        line = ''
        for u in G.neighbors(v):
            # two vertices may share multiple edges
            for data in G[u][v].values():
                line += line and ' + '
                line += 'x{}{}'.format(data['id'], data['type'] ==
                        du.ETYPE_ID and '_%s' %i or '')
        if G.nodes[v]['type'] == du.VTYPE_EXTR:
            line += ' = 2\n'
        elif G.nodes[v]['type'] == du.VTYPE_CAP:
            line += ' - 2 o{} = 0\n'.format(v)
        else:
            raise Exception('unknown node type!')
        out.write(line)


def c03(siblings, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    # siblings are specified by their unique ID
    for id1, id2 in siblings:
        out.write('x{0} - x{1} = 0\n'.format(id1, id2))


def c04(G, i, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        out.write('y{0}_{3} - y{1}_{3} + {0} x{2}{4} <= {0}\n'.format(u, v,
            data['id'], i, data['type'] == du.ETYPE_ID and '_%s' %i or ''))
        out.write('y{1}_{3} - y{0}_{3} + {1} x{2}{4} <= {1}\n'.format(u, v,
            data['id'], i, data['type'] == du.ETYPE_ID and '_%s' %i or ''))


def c05(G, i, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        if data['type'] == du.ETYPE_ID:
            out.write('y{0}_{2} + {0} x{1}_{2} <= {0}\n'.format(u, data['id'],
                i))
            out.write('y{0}_{2} + {0} x{1}_{2} <= {0}\n'.format(v, data['id'],
                i))


def c06(G, i, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for v in G.nodes():
        out.write('{0} z{0}_{1} - y{0}_{1} <= 0\n'.format(v, i))


def c07(G, i, out, genomes):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        if data['type'] == du.ETYPE_ID:
            if G.nodes[u]['id'][0] == genomes[0]:
                out.write('r{0}_{2} + x{1}_{2} <= 1\n'.format(u, data['id'],
                    i))
                out.write('r{0}_{2} + x{1}_{2} <= 1\n'.format(v, data['id'],
                    i))
            else:
                out.write('r{0}_{2} - x{1}_{2} >= 0\n'.format(u, data['id'],
                    i))
                out.write('r{0}_{2} - x{1}_{2} >= 0\n'.format(v, data['id'],
                    i))


def c08(G, i, out):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        out.write('t{2}_{3} - r{1}_{3} + r{0}_{3} - x{2}{4} >= -1\n'.format(u,
            v, data['id'], i, data['type'] == du.ETYPE_ID and '_%s' %i or ''))
        out.write('t{2}_{3} - r{1}_{3} + r{0}_{3} - x{2}{4} >= -1\n'.format(v,
            u, data['id'], i, data['type'] == du.ETYPE_ID and '_%s' %i or ''))


def c09(G, i, out, genomes):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        if data['type'] == du.ETYPE_ADJ and G.nodes[u]['id'][0] == genomes[0]:
            # do both sides of the edge:
            line = ''
            for w in G.neighbors(v):
                for data_vw in G[v][w].values():
                    if data_vw['type'] == du.ETYPE_ID:
                        line += line and ' + '
                        line += 'x{}_{}'.format(data_vw['id'], i)
            for w in G.neighbors(u):
                for data_uw in G[u][w].values():
                    if data_uw['type'] == du.ETYPE_ID:
                        line += line and ' + '
                        line += 'x{}_{}'.format(data_uw['id'], i)
            line += ' - t{}_{} >= 0\n'.format(data['id'], i)
            out.write(line)


def c10(G, i, out, genomes):
    # XXX i is not the i of the formula in the paper, but the index of the
    # pairwise comparison!
    for u, v, data in G.edges(data = True):
        if data['type'] != du.ETYPE_ADJ or G.nodes[u]['id'][0] == genomes[0]:
            out.write('t{}_{} = 0\n'.format(data['id'], i))


def c11(G, i, out):

    # by construction, telomeres are nodes with higher IDs than all other
    # nodes. We do not allow that cycles are counted using those nodes
    for v, data in G.nodes(data = True):
        if data['type'] == du.VTYPE_CAP:
            out.write('z{0}_{1} = 0\n'.format(v, i))

def c12(circ_singletons, i, out):

    for j, path in enumerate(circ_singletons.values()):
        component_vars = ['x{}{}'.format(data['id'], data['type'] ==
            du.ETYPE_ID and '_%s' %i or '') for data in path]
        print('{} - s{}_{} <= {}'.format(' + '.join(component_vars), j, i,
            len(component_vars)-1), file = out)

def c13(caps, out):

    # speed-up: let the solver know that only an even number of telomeres
    # leads to a valid solution 
    for j, (_, cap_set) in enumerate(sorted(caps.items())):
        if cap_set:
            print('{} - 2 a{} = 0'.format(' + '.join(map(lambda x: f'o{x}',
                cap_set)), j), file = out)
'''

def getAllCaps(graphs):
    res = dict((k, set()) for k in set(chain(*graphs.keys())))

    for (child, parent), G in graphs.items():

        for v, vdata in G.nodes(data=True):
            if vdata['type'] == du.VTYPE_CAP:
                res[vdata['id'][0]].add(v)
    return res


# ILP DOMAINS
#

COUNTERS = ['w','f','n','c','s','pab','pAB','pAb','pAa','pBa','pBb','pABa','pABb']
def domains(graphs, out):
    global_generals = set()
    out.write('bounds\n')
    for te, ((child, parent), G) in enumerate(sorted(graphs.items())):
        for rv in COUNTERS:
            print("0 <= {rv}{sep}{e} <= inf".format(rv=rv,sep=du.SEP,e=te),file=out)
        print("q{sep}{e} free".format(sep=du.SEP,e=te))
        for v,data in G.nodes(data=True):
            if data['type']==du.VTYPE_CAP:
                continue
            print("0 <= w{sep}{te}{sep}{v} <= inf".format(sep=du.SEP,te=te,v=v),file=out)
            global_generals.add("g{sep}{v}".format(sep=du.SEP,v=v))
    for gg in global_generals:
        print("0 <= {gg} <= inf".format(gg=gg),file=out)
        


#
# ILP VARIABLES
#

def variables(graphs, out):

    #
    # integer variables
    #
    global_generals = set()
    out.write('generals\n')
    for te, ((child, parent), G) in enumerate(sorted(graphs.items())):
        LOG.info(('writing general variables for relational diagram of {} ' + \
                'and {}').format(child, parent))
        for rv in COUNTERS:
            print("{rv}{sep}{e}".format(rv=rv,sep=du.SEP,e=te),file=out)
        print("q{sep}{e}".format(sep=du.SEP,e=te))
        for v,data in G.nodes(data=True):
            if data['type']==du.VTYPE_CAP:
                continue
            print("w{sep}{te}{sep}{v}".format(sep=du.SEP,te=te,v=v),file=out)
            print("y{sep}{te}{sep}{v}".format(sep=du.SEP,te=te,v=v),file=out)
            global_generals.add("g{sep}{v}".format(sep=du.SEP,v=v))

    for gg in global_generals:
        print(gg,file=out)
    #
    # binary variables
    #
    out.write('binaries\n')

    global_binaries = set()
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing binary variables for relational diagram of {} ' + \
                'and {}').format(child, parent))
        genomes = [child,parent]
        for v,data in G.nodes(data=True):
            print("z{sep}{v}".format(sep=du.SEP,v=v),file=out)
            print("l{sep}{v}".format(sep=du.SEP,v=v),file=out)
            if data['type']!=du.VTYPE_CAP:
                print("rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                print("d{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=v),file=out)
            else:
                if get_genome(G,v) == genomes[0]:
                    print("rAB{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rAa{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rAb{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                else:
                    print("rBa{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
        for u,v,data in G.edges(data=True):
            print("x{sep}{te}{sep}{e}".format(sep=du.SEP,te=tree_edge,e=data['id']),file=out)

    print('\n'.join(global_binaries), file = out)
    print('\n')


def identifyCandidateTelomeres(candidateAdjacencies, weightThreshold, dont_add=False):

    res = dict()
    weights = candidateAdjacencies['weights']
    for species, adjs in candidateAdjacencies['adjacencies'].items():
        genes = candidateAdjacencies['genes'][species]
        # add gene extremities incident to telomeric adjacencies to telomere
        # set
        telomeres = set(x[0][0] for x in adjs if x[0][1] == 'o').union(
                (x[1][0] for x in adjs if x[1][1] == 'o'))

        # remove telomeric extremities from gene set
        for t in telomeres:
            genes.remove(t)

        if not dont_add:
            G = nx.Graph()
            G.add_nodes_from(reduce(lambda x, y: x + y, (((g, du.EXTR_HEAD), (g,
                du.EXTR_TAIL)) for g in genes)))
            G.add_edges_from(adjs)

            for C in tuple(nx.connected_components(G)):
                C = set(C)

                # check if component is linear / circular / fully connected /
                # even odd
                # - if it is (linear/circular or fully connected) and even, no
                # telomere needs to be added
                degs = set(map(lambda x: x[1], G.degree(C)))

                # structures that support a perfect matching, e.g.
                # - simple paths/simple cycles of even size
                # - fully connected components of even size
                # do not need telomeres and are omitted by the condition below.
                # Conversely, nodes of components of odd size (including components
                # of size 1) will always be considered as candidate telomeres

                # will evaluate to true if component is NOT
                # - linear/circular or
                # - fully connected
                # - even (in terms of #nodes)
                if degs.difference((1, 2)) and degs.difference((len(C)-1,)) or len(C) % 2:
                    for g, extr in C:
                        t = f't_{g}_{extr}'
                        telomeres.add(t)
                        adjs.append(((g, extr), (t, 'o')))

#        genes_edg = [((g, du.EXTR_HEAD), (g, du.EXTR_TAIL)) for g in genes]
#        if species == 'n3':
#            C = nx.connected.node_connected_component(G, ('69_6', 'h'))
#            G = G.subgraph(C).copy()
##            G.add_edges_from(genes_edg)
#            pos = nx.spring_layout(G)
#            #nx.draw_networkx_nodes(G, pos=pos, node_size=8)
#            nx.draw(G, pos=pos, node_size=10)
#            #nx.draw_networkx_labels(G, pos=pos, font_size=12)
#            nx.draw_networkx_edges(G, pos, set(map(tuple,
#                adjs)).intersection(G.edges()), edge_color='black')
##            nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=dict(((x[0], \
##                    x[1]), G[x[0]][x[1]][0]['type']) for x in G.edges(data = \
##                    True)))
##            nx.draw_networkx_edges(G, pos, set(genes_edg).intersection(G.edges()),
##                edge_color='red')
#            import matplotlib.pylab as plt
#            import pdb; pdb.set_trace()
        res[species] = telomeres
        LOG.info('identified %s candidate telomeres in genome %s' %(
            len(telomeres), species))
    return res


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=open,
            help='phylogenetic tree as parent-child relation table')
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('-t', '--no_telomeres', action='store_true',
            help='don\'t add any additional telomeres')
    parser.add_argument('-m', '--output_id_mapping', type=FileType('w'),
            help='writs a table with ID-to-gene extremity mapping')
    parser.add_argument('-a', '--alpha', default = 0.5, type=float,
            help='linear weighting factor for adjacency weights vs DCJ ' + \
                    'indel distance (alpha = 1 => maximize only DCJ indel ' + \
                    'distance)')
    parser.add_argument('-b', '--beta', default = -1, type=float,
            help='linear weighting factor for telomeric adjacencies;' + \
                    'if beta < 0, then beta is set to 1/2 * alpha')
    parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    beta = args.beta
    if beta < 0:
        beta = args.alpha * 0.5

    # load & process input data
    LOG.info('loading species tree from {}'.format(args.tree.name))
    speciesTree = du.parseTree(args.tree)

    LOG.info(('loading candidate adjacencies from {}, using "{}" to separate' + \
            ' gene family from uniquifying identifier').format(
        args.candidateAdjacencies.name, args.separator))
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)

    # add telomeres
    telomeres = identifyCandidateTelomeres(candidateAdjacencies,
            ADJ_TRUST_THRESHOLD, args.no_telomeres)

    # construct adjacency graphs
    genes = candidateAdjacencies['genes']
    adjacencies = candidateAdjacencies['adjacencies']
    weights = candidateAdjacencies['weights']
    penalities = candidateAdjacencies['penalities']

    ext2id = du.IdManager()
    LOG.info(('constructing relational diagrams for all {} branches of ' + \
            'the tree').format(len(speciesTree)))
    relationalDiagrams = du.constructRelationalDiagrams(speciesTree,
            adjacencies, telomeres, weights, penalities, genes, ext2id,
            sep=args.separator)

    graphs = relationalDiagrams['graphs']

#    for gNames, G in graphs.items():
#            genes_edg = list()
#            for gName in gNames:
#                genes = candidateAdjacencies['genes'][gName]
#                genes_edg.extend(((ext2id.getId((gName, (g, du.EXTR_HEAD))),
#                    ext2id.getId((gName, (g, du.EXTR_TAIL))))for g in genes))
##            Gp = nx.Graph()
##            Gp.add_edges_from(genes_edg)
##            Gp.add_edges_from((u, v) for u, v, data in G.edges(data=True) if
##                    data['type'] == du.ETYPE_ADJ)
#            G = G.subgraph(nx.node_connected_component(G, ext2id.getId(('n10',
#                ('23_7', 'h')))))
#            pos = nx.spring_layout(G)
##            nx.draw_networkx_edges(G, pos, set(genes_edg).intersection(G.edges()),
##                edge_color='red')
#            nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in
#                G.edges(data=True) if data['type'] == du.ETYPE_EXTR],
#                edge_color='green')
#            nx.draw_networkx_nodes(G, pos=pos, node_size=8)
#            #nx.draw(G, pos=pos, node_size=10)
#            nx.draw_networkx_labels(G, pos=pos, font_size=10, labels = dict((v,
#                '{0}:{1[0]}{1[1]}'.format(*G.nodes[v]['id'])) for v in G.nodes()))
#            nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in
#                G.edges(data=True) if data['type'] == du.ETYPE_ID],
#                edge_color='gray')
#            nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in
#                G.edges(data=True) if data['type'] == du.ETYPE_ADJ],
#                edge_color='red')
##            nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=dict(((x[0], \
##                    x[1]), G[x[0]][x[1]][0]['type']) for x in G.edges(data = \
##                    True)))
#            import matplotlib.pylab as plt
#            import pdb; pdb.set_trace()

    siblings = relationalDiagrams['siblings']

    for G in graphs.values():
        du.checkGraph(G)

    #circ_singletons = dict()
    #for ident, G in graphs.items():
    #    circ_singletons[ident] = du.identifyCircularSingletonCandidates(G)
    #    LOG.info(f'identified {len(circ_singletons[ident])} circular singleton candidates')

    #caps = getAllCaps(graphs)
    # construct & output ILP
    out = stdout

    LOG.info('writing objective over all graphs')
    objective(graphs, args.alpha, out)

    LOG.info('writing constraints...')
    constraints(graphs, siblings, out)

    LOG.info('writing domains...')
    domains(graphs, out)

    LOG.info('writing variables...')
    variables(graphs, out)

    if args.output_id_mapping:
        LOG.info('writing ID-to-gene extremity mapping to {}'.format(
            args.output_id_mapping.name))
        idMap = ext2id.getMap()
        out_table = list()
        for k, v in idMap.items():
            out_table.append((str(v), k[0], k[1][0], k[1][1]))
        out_table.sort(key = lambda x: int(x[0]))
        print('\n'.join(map(lambda x: '\t'.join(x), out_table)),
                file=args.output_id_mapping)

    LOG.info('DONE')
    out.write('end\n')


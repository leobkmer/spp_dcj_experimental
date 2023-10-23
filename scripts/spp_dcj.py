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

def objective(graphs, circ_singletons, alpha, beta, out):
    out.write('maximize ')

    # sum of adjacency weights over all genomes
    adjs = set(reduce(lambda x, y: x + y, (tuple(map(lambda z: (z[2]['id'], \
            z[2]['weight']), filter(lambda x: x[2]['type'] == du.ETYPE_ADJ, \
            G.edges(data = True)))) for i, (_, G) in \
            enumerate(sorted(graphs.items())))))

    out.write(' + '.join(map(lambda x: '{} x{}'.format(x[1] * (1-alpha), \
            x[0]), adjs)))

    for i, (_, G) in enumerate(sorted(graphs.items())):
        if G.number_of_edges():
            out.write(' + ')
        # DCJ indel distance
        out.write(' + '.join(map(lambda x: '{} z{}_{}'.format(alpha, x, i),
            G.nodes())))
        out.write(''.join(map(lambda x: ' - {} t{}_{}'.format(0.5 * alpha,
            x[2]['id'], i), G.edges(data = True))))

    for i, ident in enumerate(sorted(graphs.keys())):
        cs = circ_singletons[ident]
        if cs:
            out.write(' - ')
        # subtract circular singleton penality
        out.write(' - '.join(('{} s{}_{}'.format(alpha, j, i) for j in range(len(cs)))))

    for adj, weight in set(reduce(lambda x, y: x + y, (tuple(map(lambda z: (z[2]['id'], \
            z[2]['penality']), filter(lambda x: 'penality' in x[2] and \
            x[2]['type'] == du.ETYPE_ADJ, G.edges(data = True)))) for i, (_, G) in \
            enumerate(sorted(graphs.items()))))):
        out.write(f' - {weight * beta} x{adj}')

    out.write('\n\n')


#
# ILP CONSTRAINTS
#

def constraints(graphs, siblings, circ_singletons, caps, out):

    out.write('subject to\n')

    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing constraints for relational diagram of {} and ' + \
                '{}').format(child, parent))
        c01(G, out)
        c02(G, i, out)
        c03(siblings[(child, parent)], out)
        c04(G, i, out)
        c05(G, i, out)
        c06(G, i, out)
        c07(G, i, out, (child, parent))
        c08(G, i, out)
        c09(G, i, out, (child, parent))
        c10(G, i, out, (child, parent))
        c11(G, i, out)
        c12(circ_singletons[(child, parent)], i, out)
        c13(caps, out)

    out.write('\n')





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




DEFAULT_SEP = '_'
ENUMERATE_CONSTRAINTS = False

#give a canonical string for a given edge between u,v with index k
def disp_edge(u,v,k,sep=DEFAULT_SEP):
    p = sep.join(sorted([u,v]))
    return '{p}{s}{k}'.format(p=p,s=sep,k=k)



'''
I presume here that pseudo-caps in the graph are declared as du.VTYPE_CAP
'''

# implement constraints C.04 to C.37
def mrd_constraints(graphs, siblings, circ_singletons, pseudocaps, out):

    out.write('subject to\n')

    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing constraints for relational diagram of {} and ' + \
                '{}').format(child, parent))
        
        
        

    out.write('\n')


def cfc04(i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    print("{enm}f{s}{e} - n{s}{e} + c{s}{e} - q{s}{e} - s{s}{e} = 0".format(s=sep,e=i,enm="c04: " if enm else ""),file=out)

def cfc05(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['0.5 x{s}{i}{s}{e}'.format(s=sep,e=disp_edge(u,v,k),i=i) for u,v,k,tp in G.edges(keys=True,data=type) if tp==du.ETYPE_EXTR])
    print("{enm}{sumstr} - n{s}{e} = 0".format(sumstr=sumstr,e=i,s=sep,enm="c05: " if enm else ""),file=out)

def cfc06(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rc{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,vt in G.nodes(data='type') if vt != du.VTYPE_CAP])
    print("{enm}{sumstr} - c{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm= "c.06: " if enm else ""),file=out)

def cfc07(i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    print("{enm}piiab{s}{i} + pABa{s}{i} + pABb{s}{i} - pttAB{s}{i} - 2 q{s}{i} <= 0".format(i=i,s=sep,enm="c.07: " if enm else ""),file=out)

def cfc08(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['riiab{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP and data['genome'] == du.VGENOME_FIRST])
    print("{enm}{sumstr} - piiab{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.08: " if enm else ""),file=out)

def cfc09(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rtiAb{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] == du.VTYPE_CAP and data['genome'] == du.VGENOME_FIRST])
    print("{enm}{sumstr} - ptiAb{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.09: " if enm else ""),file=out)

def cfc10(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rtiBa{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] == du.VTYPE_CAP and data['genome'] == du.VGENOME_SECOND])
    print("{enm}{sumstr} - ptiBa{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.09: " if enm else ""),file=out)


def cfc11(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rtiAa{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] == du.VTYPE_CAP and data['genome'] == du.VGENOME_FIRST])
    print("{enm}{sumstr} - ptiAa{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.09: " if enm else ""),file=out)

def cfc12(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rtiBb{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] == du.VTYPE_CAP and data['genome'] == du.VGENOME_SECOND])
    print("{enm}{sumstr} - ptiBb{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.09: " if enm else ""),file=out)


def cfc13to16(i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    count = 13
    for x in 'ab':
        for y in 'BA':
            print("{enm}pAB{x}{s}{i} - pti{y}{x} >= 0".format(x=x,y=y,s=sep,i=i,enm="c.%02d: "%count if enm else ""),file=out)
            count+=1

def cfc17(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rttAB{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] == du.VTYPE_CAP and data['genome'] == du.VGENOME_FIRST])
    print("{enm}{sumstr} - pttAB{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.09: " if enm else ""),file=out)


def cfc18(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    sumstr = ' + '.join(['rs{s}{i}{s}{v}'.format(i=i,v=v,s=sep) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP])
    print("{enm}{sumstr} - s{s}{i} = 0".format(sumstr=sumstr,i=i,s=sep,enm="c.08: " if enm else ""),file=out)

def cfc19(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc, v in enumerate(G.nodes(),start=1):
        sumstr = ' + '.join(['x{s}{i}{s}{e}'.format(i=i,s=sep,e=disp_edge(u,v,k)) for u in G[v] for k in G[v][u] if G[v][u][k] in [du.ETYPE_ADJ,du.ETYPE_EXTR]])
        print("{enm}{sumstr} - g{s}{v} = 0".format(sumstr,v=v,s=sep,enm="c.19.%d: "%cc if enm else ""),file=out)


def cfc20(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc, e in enumerate([disp_edge(u,v,k) for u,v,k,tp in G.edges(keys=True,data='type') if tp==du.ETYPE_ADJ],start=1):
        print("{enm}a{s}{e} - x{s}{i}{s}{e}".format(e=e,i=i,s=sep,enm="c.20.%d: "%cc if enm else ""),file=out)


def cfc21(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc, v in enumerate(G.nodes(),start=1):
        print("{enm}z{s}{i}{s}{v} - g{s}{v} <= 0".format(i=i,v=v,s=sep,enm="c.21.%d: "%cc if enm else ""),file=out)


def cfc22(sibs,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc,es in enumerate(sibs):
        e1 = disp_edge(*es[0])
        e2 = disp_edge(*es[1])
        print("{enm}x{s}{i}{s}{e1} - x{s}{i}{s}{e2} = 0".format(i=i,e1=e1,e2=e2,s=sep,enm="c22.%d: "%cc if enm else ""),file=out)

#TODO: implement this after checking w/ dany
def cfc23(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for pacde in G.edges(keys=True,data=True):
        u_,v_,k,data = pacde
        e = disp_edge(u_,v_,k)
        for u,v in [(u_,v_),(v_,u_)]:
            if data['type'] != du.ETYPE_ID:
                print("{enm}y{s}{i}{s}{u} - y{s}{i}{s}{v} + {ixu} x{s}{i}{s}{e}  <= {ixu}".format(
                    i=i,u=u,v=v,e=e,ixu=G.nodes[u]['localid'],s=sep,enm='c23.%d '%cc if enm else ""
                ),file=out)
                cc+=1
            else:
                print("{enm}y{s}{i}{s}{u} + {ixu} x{s}{i}{s}{e}  <= {ixu}".format(
                    i=i,u=u,e=e,ixu=G.nodes[u]['localid'],s=sep,enm='c23.%d '%cc if enm else ""
                ),file=out)
                cc+=1

def cfc24(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc,pacdv in enumerate(G.nodes(data='localid')):
        v,lid = pacdv
        print("{enm}{lid} z{s}{i}{s}{v} - y{s}{i}{s}{v} <= 0".format(i=i,v=v,lid=lid,s=sep,enm="c24.%d: "%cc if enm else ''),file=out)

def cfc25(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for pacde in [(u,v,k) for u,v,k,et in G.edges(keys=True,data='type') if et==du.ETYPE_ID]:
        u_,v_,k = pacde
        e = disp_edge(u_,v_,k)
        for v in [u_,v_]:
            if G.nodes[v]['genome'] == du.VGENOME_FIRST:
                print("{enm}l{s}{i}{s}{v} + x{s}{i}{s}{e} <= 1".format(i=i,s=sep,v=v,e=e,enm="c25.%d: "%cc if enm else ""),file=out)
            else:
                print("{enm}x{s}{i}{s}{e} - l{s}{i}{s}{v} <= 0".format(i=i,s=sep,v=v,e=e,enm="c25.%d: "%cc if enm else ""),file=out)
            cc+=1

def cfc26(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for u_,v_,k,et in G.edges(keys=True,data='type'):
        if et not in [du.ETYPE_ADJ,du.ETYPE_EXTR]: #and du.VTYPE_CAP not in [G.nodes[u_]['type'],G.nodes[v_]['type']]:
            continue
        e = disp_edge(u_,v_,k)
        for u,v in[u_,v_]:
            if et == du.ETYPE_EXTR:
                print("{enm}l{s}{i}{s}{v} - l{s}{i}{s}{u} + x{s}{i}{s}{e} <= 1".format(i=i,v=v,u=u,e=e,s=sep,enm="c25.%d: "%cc if enm else ""),
                    file=out)
                cc+=1
            elif et == du.ETYPE_ADJ:
                if du.VTYPE_CAP in [G.nodes[u_]['type'],G.nodes[v_]['type']]:
                    continue
                if G.nodes[v]['genome']==du.VGENOME_FIRST:
                    print("{enm}l{s}{i}{s}{v} - l{s}{i}{s}{u} + x{s}{i}{s}{e} - riiab{s}{i}{s}{v} - riiab{s}{i}{s}{u} <= 1".format(i=i,v=v,u=u,e=e,s=sep,enm="c25.%d: "%cc if enm else ""),
                    file=out)
                else:
                    print("{enm}l{s}{i}{s}{v} - l{s}{i}{s}{u} + x{s}{i}{s}{e} <= 1".format(i=i,v=v,u=u,e=e,s=sep,enm="c26.%d: "%cc if enm else ""),
                    file=out)
                cc+=1

def cfc27(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for cc,v in enumerate([v for v,tp in G.nodes(data='type') if tp != du.VTYPE_CAP]):
        print("{enm}rc{s}{i}{s}{v} - z{s}{i}{s}{v} <= 0".format(i=i,v=v,s=sep,enm="c27.%d: "%cc if enm else ""),file=out)

def cfc28(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for u_,v_,k,tp in G.nodes(keys=True,data='type'):
        if tp != du.ETYPE_ID:
            continue
        e = disp_edge(u_,v_,k)
        for v in[u_,v_]:
            print("{enm}riiab{s}{i}{s}{v} - x{s}{i}{s}{e} <= 0".format(i=i,v=v,e=e,s=sep,enm="c28.%d: "%cc if enm else ""),file=out)
        cc+=1

def cfc29(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        g = 0 if data['genome'] == du.VGENOME_FIRST else 1
        print("{enm}l{s}{i}{s}{v} = {g}".format(i=i,v=v,s=sep,enm="c29.%d: "%cc if enm else ""),file=out)

def cfc30(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc = 1
    for v_,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        adje = [(v_,u,k) for u in G[v_] for k in G[v_][u] if G[v_][u][k]['type'] == du.ETYPE_ADJ]
        if len(adje) != 1:
            raise AssertionError("Pseudo-cap connected to %d vertices."%len(adje))
        if data['genome']==du.VGENOME_FIRST:
            v,u,k = adje[0]
            e = disp_edge(v,u,k)
            print("{enm}l{s}{i}{s}{u} - l{s}{i}{s}{v} - rttAB{s}{i}{s}{v} - rtiAb{s}{i}{s}{v} + x{s}{i}{s}{e} <= 1".format(
                i=i,v=v,u=u,e=e,s=sep,enm="c30.%d: "%cc if enm else ""),file=out)
            cc+=1
        else:
            u,v,k = adje[0]
            e = disp_edge(v,u,k)
            print("{enm}l{s}{i}{s}{u} - l{s}{i}{s}{v} - rtiBa{s}{i}{s}{v} + x{s}{i}{s}{e} <= 1".format(
                i=i,v=v,u=u,e=e,s=sep,enm="c30.%d: "%cc if enm else ""),file=out)
            cc+=1
        

def cfc31(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc=1
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP or data['genome']!= du.VGENOME_FIRST:
            continue
        print("{enm}rttAB{s}{i}{s}{v} - z{s}{i}{s}{v} <= 0".format(i=i,v=v,s=sep,enm="c31.%d: "%cc if enm else ""),file=out)
        cc+=1

def cfc32(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc=1
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        rvars = ['rtiAb','rtiAa'] if data['genome']==du.VGENOME_FIRST else ['rtiBa','rtiBb']
        print("{enm}{ra}{s}{i}{s}{v} + {rb}{s}{i}{s}{v} + y{s}{i}{s}{v} >= 1".format(
            ra=rvars[0],rb=rvars[1],v=v,i=i,s=sep,enm="c32.%d: "%cc if enm else ""),file=out)
        cc+=1

def cfc33(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        rvars = ['rtiAb','rtiAa'] if data['genome']==du.VGENOME_FIRST else ['rtiBa','rtiBb']
        for r in rvars:
            print("{enm}y{s}{i}{s}{v} + {ixv} {r}{s}{i}{s}{v} <= {ixv}".format(i=i,v=v,s=sep,ixv=data['localid'],enm="c33.%d: "%cc if enm else ""),file=out)
            cc+=1

def cfc34(G,i,out,sep=DEFAULT_SEP,enm=ENUMERATE_CONSTRAINTS):
    cc=1
    for v_,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        adje = [(v_,u,k) for u in G[v_] for k in G[v_][u] if G[v_][u][k]['type'] == du.ETYPE_ADJ]
        if len(adje) != 1:
            raise AssertionError("Pseudo-cap connected to %d vertices."%len(adje))
        v,u,k = adje[0]
        if data['genome']==du.VGENOME_FIRST:
            print("{enm}rttAB{s}{i}{s}{v} - l{s}{i}{s}{u} <= 0".format(i=i,v=v,u=u,s=sep,enm="c34.%d: "%cc if enm else ""))
            cc+=1
            print("{enm}rtiAb{s}{i}{s}{v} - l{s}{i}{s}{u} <= 0".format(i=i,v=v,u=u,s=sep,enm="c34.%d: "%cc if enm else ""))
            cc+=1
        else:
            print("{enm}rtiBa{s}{i}{s}{v} + l{s}{i}{s}{u} <= 1".format(i=i,v=v,u=u,s=sep,enm="c34.%d: "%cc if enm else ""))
            cc+=1

def getAllCaps(graphs):
    res = dict((k, set()) for k in set(chain(*graphs.keys())))

    for (child, parent), G in graphs.items():

        for v, vdata in G.nodes(data=True):
            if vdata['type'] == du.VTYPE_CAP:
                res[vdata['id'][0]].add(v)
    return res


# ILP DOMAINS
#

def domains(graphs, out):

    out.write('bounds\n')

    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):
        LOG.info(('writing domains for relational diagram of {} and ' + \
                '{}').format(child, parent))
        d02(G, i, out)

    out.write('\n')


def d02(G, i, out):

    for v in G.nodes():
        out.write('0 <= y{0}_{1} <= {0}\n'.format(v, i))


#
# ILP VARIABLES
#

def variables(graphs, circ_singletons, caps, out):

    #
    # integer variables
    #
    out.write('generals\n')

    variables = set()
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing general variables for relational diagram of {} ' + \
                'and {}').format(child, parent))

        # D.02
        for v in G.nodes():
            variables.add('y{}_{}'.format(v, i))


    # D.08
    for j, (_, cap_set) in enumerate(sorted(caps.items())):
        if cap_set:
            variables.add(f'a{j}')

    print('\n'.join(variables), file = out)
    print('\n')
    #
    # binary variables
    #
    out.write('binaries\n')

    variables = set()
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing binary variables for relational diagram of {} ' + \
                'and {}').format(child, parent))

        # D.01
        for _, _, data in G.edges(data = True):
            variables.add('x{}{}'.format(data['id'], data['type'] ==
                du.ETYPE_ID and '_%s' %i or ''))

        # D.03
        for v in G.nodes():
            variables.add('z{}_{}'.format(v, i))

        # D.04
        for v in G.nodes():
            variables.add('r{}_{}'.format(v, i))

        # D.05
        for _, _, data in G.edges(data = True):
            variables.add('t{}_{}'.format(data['id'], i))

        # D.06
        for v, type_ in G.nodes(data='type'):
            if type_ == du.VTYPE_CAP:
                variables.add('o{}'.format(v))

        # D.07
        for j in range(len(circ_singletons[(child, parent)])):
            variables.add('s{}_{}'.format(j, i))

    print('\n'.join(variables), file = out)
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

    LOG.info('loading candidate adjacencies from {}'.format(
        args.candidateAdjacencies.name))
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies)

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
            adjacencies, telomeres, weights, penalities, genes, ext2id)

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

    circ_singletons = dict()
    for ident, G in graphs.items():
        circ_singletons[ident] = du.identifyCircularSingletonCandidates(G)
        LOG.info(f'identified {len(circ_singletons[ident])} circular singleton candidates')

    caps = getAllCaps(graphs)
    # construct & output ILP
    out = stdout

    LOG.info('writing objective over all graphs')
    objective(graphs, circ_singletons, args.alpha, beta, out)

    LOG.info('writing constraints...')
    constraints(graphs, siblings, circ_singletons, caps, out)

    LOG.info('writing domains...')
    domains(graphs, out)

    LOG.info('writing variables...')
    variables(graphs, circ_singletons, caps, out)

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


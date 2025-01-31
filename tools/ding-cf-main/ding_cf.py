#!/usr/bin/python3
from argparse import ArgumentParser, FileType, ArgumentTypeError
import argparse as ap
import logging
from ilp_util_adj import *
import networkx as nx
import math
import sys
from dingII_util import *
from copy import deepcopy

MAX_CONSTRAINTS = 100
SHOW_LEGACY_VERSIONS = False

def matching_range(v):
    iv = float(v)
    if iv < 0.0 or iv > 1.0:
        raise ArgumentTypeError("%s is an invalid range bound. Please only provide values in [0,1]." % v)
    return iv

def mm_bound(a,b):
    m = min(a,b)
    return (m,m)
    
def em_bound(a,b):
    if a == 0 or b == 0:
        return (0,0)
    return (1,1)

def im_bound(a,b):
    if a == 0 or b == 0:
        return (0,0)
    m = min(a,b)
    return (1,m)

class RangedBound:
    def r_bound(self, a, b):
        m = min(a,b)
        return (int(math.ceil(self.lower*m)), int(math.ceil(self.upper*m)))



def create_model(amults, bmults, boundf):
    bounds = {}
    for gene, mult in amults.items():
        bounds[gene] = boundf(mult, dgetc(bmults, gene))
    for gene, mult in bmults.items():
        bounds[gene] = boundf(mult, dgetc(amults, gene))
    return bounds
        

def create_matching_model(args, gnms):
    amults = get_multiplicities(gnms[0])
    bmults = get_multiplicities(gnms[1])
    needs_supplement = True
    model = {}
    if args.custom:
        model = read_matching_model(args.custom)
        enforce_model_integrity(model, amults, bmults)
        needs_supplement = False
        if is_incomplete(model, amults, bmults):
            LOG.warning('The matching model provided via custom file is incomplete. It will be supplemented by your chosen model.')
            needs_supplement = True
    if needs_supplement:    
        if args.exemplary:
            boundf = em_bound
        elif args.intermediate:
            boundf = im_bound
        elif args.range:
            rb = RangedBound()
            bs = sorted(args.range)
            rb.lower = bs[0]
            rb.upper = bs[1]
            boundf = rb.r_bound
        else:
            if not args.maximal:
                LOG.info('No matching model chosen. Default is maximal matching.')
            boundf = mm_bound
        model2 = create_model(amults, bmults, boundf)
        #supplement model so far
        for gene in model2:
            if not gene in model:
                model[gene] = model2[gene]
    return model
    
def canon_disp_tup(tp):
    tp_ = tp[0:2]
    return '%i_%i_%i'%(min(tp_),max(tp_),tp[2])

def get_multiplicities(gnm):
    mlt = {}
    for _, chrm in gnm[1]:
        for _, gene in chrm:
            dadd(mlt,gene,1)
    return mlt


def write_matching_model(model, fl):
    fl.write('Gene\tLower\tUpper\n')
    for gene, bnds in model.items():
        l, u = bnds
        fl.write('%s\t%i\t%i\n'%(gene,l,u))

def read_matching_model(fl):
    lines = fl.readlines()
    #throw away header
    lines = lines[1:]
    model = {}
    for line in lines:
        entries = line.strip().split('\t')
        model[entries[0]] = (int(entries[1]),int(entries[2]))
    return model

def enforce_model_integrity(model, amults, bmults):
    for gene, bnds in model.items():
        l, u = bnds
        m = min(dgetc(amults, gene),dgetc(bmults, gene))
        if l < 0:
            LOG.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if u < 0:
            LOG.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if l > u:
            LOG.warning('Inconsistent bounds: upper below lower (lower %i, upper %i) for gene %s. Both bounds will be set to %i.'%(l,u,gene, u))
            l = u
        if l > m:
            LOG.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(l,gene,m,m))
            l = m
        if u > m:
            LOG.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(u,gene,m,m))
            u = m
        if l == 0 and m > 0:
            LOG.info('Lower bound found for gene %s is 0. If this affects many genes and is not counteracted by deletion penalties this may lead to very large indel stretches.'%gene)
        model[gene] = (l,u)

def is_incomplete(model, amults, bmults):
    for gene in amults:
        if not gene in model:
            return True
    for gene in bmults:
        if not gene in model:
            return True
    return False

def get_first_extremity(orient):
    if orient == ORIENT_POSITIVE:
        return TAIL
    else:
        return HEAD

def get_last_extremity(orient):
    if orient == ORIENT_POSITIVE:
        return HEAD
    else:
        return TAIL



def siblings(rd, gi1, gi2):
    sibs = {}
    for g, occs1 in gi1[HEAD].items():
        occs2 = getl(gi2[HEAD],g)
        for hid1 in occs1:
            for hid2 in occs2:
                edgeid1 = list(rd[hid1][hid2].keys())
                if len(edgeid1) != 1:
                    raise Exception('Wrong number of extremity edges (%i) between vertices %i %i.'%(len(edgeid1),hid1,hid2))
                sister = (hid1, hid2,edgeid1[0])
                tid1 = maybe_adjacent_selfedge_raw(rd,hid1)
                if len(tid1) != 1:
                    raise Exception('Vertex %i connected to %i self edges!'%(hid1, len(tid1)))
                else:
                    tid1 = tid1[0][1]
                tid2 = maybe_adjacent_selfedge_raw(rd,hid2)
                if len(tid2) != 1:
                    raise Exception('Vertex %i connected to %i self edges!'%(hid2, len(tid2)))
                else:
                    tid2 = tid2[0][1]
                edgeid2 = list(rd[tid1][tid2].keys())
                if len(edgeid2) != 1:
                    raise Exception('Wrong number of extremity edges (%i) between vertices %i %i.'%(len(edgeid2),tid1,tid2))
                brother = (tid1, tid2,edgeid2[0])
                insertl(sibs, g, (brother, sister))
    return sibs
                

def add_constraint(ilp, constraint_type, constraint_string):
    cid = len(getl(ilp, constraint_type))
    insertl(ilp, constraint_type, (cid, constraint_string))
    

def maybe_adjacent_selfedge(rd, u):
    return [canon_disp_tup(t) for t in maybe_adjacent_selfedge_raw(rd, u)]
    
def maybe_adjacent_selfedge_raw(rd,u):
    return [(u,v,k) for v, edges in rd[u].items() for k, data in edges.items() if data['etype']==SELFEDGE]

def edge_constraints(rd, ilp,bm=False):
    '''
    Add constraint by iterating through edges once.
    '''
    for u,v, k, data in rd.edges(data=True, keys=True):
        e = canon_disp_tup((u,v,k))
        etype= data['etype']
        if etype == SELFEDGE:
            insertl(ilp, 'binary', 'd_%s'%e)
            if rd.nodes[u]['genome'] == GENOME1:
                icode = code['a']
            elif rd.nodes[u]['genome'] == GENOME2:
                icode = code['b']
            for i in [u,v]:
                if not bm:
                    #c.08
                    cs = 'm_{i} + {eightminusK} d_{e} <= 8'.format(i=i,eightminusK=8-icode,e=e)
                    add_constraint(ilp,'c08',cs)
                    #c.09
                    cs = 'm_{i} - {K} d_{e} >= 0'.format(i=i,K=icode,e=e)
                    add_constraint(ilp,'c09',cs)
                else:
                    #c.08
                    x = 'a' if rd.nodes[u]['genome'] == GENOME1 else 'b'
                    cs = 'm_{i}_{x} - d_{e} >= 0'.format(i=i,e=e,x=x)
                    add_constraint(ilp,'c08',cs)
            continue
        insertl(ilp, 'binary', 'x_%s'%e)
        if etype == ADJACENCY:
            #c.01
            cs = 'x_%s = 1'%e
            add_constraint(ilp, 'c01', cs)
        #c.10
        if etype==ADJACENCY:
            messengerv = 'n'
        else:
            messengerv = 'm'
        for i,j in [(u,v),(v,u)]:
            if not bm:
                cs = '{mess}_{i} - {mess}_{j} + 8 x_{e} <= 8'.format(mess=messengerv,i=i,j=j,e=e)
                add_constraint(ilp,'c10',cs)
            else:
                 for x in SINGLE_LABS:
                     cs =  '{mess}_{i}_{x} - {mess}_{j}_{x} + x_{e} <= 1'.format(mess=messengerv,i=i,j=j,e=e,x=x)
                     add_constraint(ilp,'c10',cs)
        #c.05
        for i,j in [(u,v), (v,u)]:
            cs = 'y_%i - y_%i + %i x_%s <= %i'%(i,j,i,e, i)
            add_constraint(ilp, 'c05', cs)
        if etype == EXTREMITYEDGE:
            insertl(ilp, 'obj', (0.5, 'x_%s'%e))


R_AB = 'ttAB'
R_Ab = 'tiAb'
R_Aa = 'tiAa'
R_ab = 'iiab'
R_Ba = 'tiBa'
R_Bb = 'tiBb'
R_C = 'C'

def any_caps(l,rd):
    return TELO in [rd.nodes[v]['extremity'] for v in l]

def edge_constraints_n(rd, ilp):
    '''
    Add constraint by iterating through edges once.
    '''
    for u,v, k, data in rd.edges(data=True, keys=True):
        e = canon_disp_tup((u,v,k))
        etype= data['etype']
        #add edge variables
        insertl(ilp,'binary','x_{e}'.format(e=e))
        if etype==ADJACENCY:
            if not any_caps([u,v],rd):
                insertl(ilp,'binary','r{c}_{e}'.format(c=R_C,e=e))
                if rd.nodes[v]==GENOME1:
                    insertl(ilp,'binary','r{ab}_{e}'.format(e=e,ab=R_ab))
            else:
                if rd.nodes[v]['genome'] == GENOME1:
                    for r in [R_Aa,R_Ab,R_AB]:
                        insertl(ilp,'binary','r{r}_{e}'.format(r=r,e=e))
                else:
                    for r in [R_Ba,R_Bb]:
                        insertl(ilp,'binary','r{r}_{e}'.format(r=r,e=e))
        if etype==ADJACENCY:
            add_constraint(ilp,'c01',"x_{e} = 1".format(e=e))
        for i,j in [(u,v),(v,u)]:
            if etype in [ADJACENCY,EXTREMITYEDGE]:
                add_constraint(ilp,'c04',"y_{j} - y_{i} + {j} x_{e} <= {j}".format(i=i,j=j,e=e))
            elif etype == SELFEDGE:
                add_constraint(ilp,'c04',"y_{j} + {j} x_{e} <= {j}".format(j=j,e=e))
            if etype == SELFEDGE:
                if rd.nodes[i]['genome'] == GENOME1:
                    add_constraint(ilp,'c08',"b_{i} + x_{e} <= 1".format(i=i,e=e))
                else:
                    add_constraint(ilp,'c08',"b_{i} - x_{e} >= 0".format(i=i,e=e))
            if etype==EXTREMITYEDGE:
                add_constraint(ilp,'c09',"b_{i} - b_{j} + x_{e} <= 1".format(i=i,j=j,e=e))
            elif etype == ADJACENCY:
                if not any_caps([i,j],rd) and rd.nodes[v]['genome']==GENOME1:
                    add_constraint(ilp,'c09',"b_{i} - b_{j} - r{ab}_{e} + x_{e} <= 1".format(i=i,j=j,ab=R_ab,e=e))
                elif rd.nodes[j]['genome'] == GENOME1 and rd.nodes[j]['extremity']==TELO:
                    add_constraint(ilp,'c09',"b_{i} - b_{j} + x_{e} - r{Ab}_{e} - r{AB}_{e}  <= 1".format(i=i,j=j,e=e,Ab=R_Ab,AB=R_AB))
                elif rd.nodes[i]['genome'] == GENOME2 and rd.nodes[i]['extremity']==TELO:
                    add_constraint(ilp,'c09',"b_{i} - b_{j} + x_{e} - r{Ba}_{e} <= 1".format(i=i,j=j,e=e,Ba=R_Ba))
                else:
                    add_constraint(ilp,'c09',"b_{i} - b_{j} + x_{e} <= 1".format(i=i,j=j,e=e))
        if etype == ADJACENCY:
            if not any_caps([u,v],rd):
                    add_constraint(ilp,'c10',"r{c}_{e} - z_{u} - z_{v} <= 0".format(u=u,v=v,e=e,c=R_C))
            elif rd.nodes[v]['genome']==GENOME1:
                    add_constraint(ilp,'c10',"r{AB}_{e} - z_{u} - z_{v} <= 0".format(u=u,v=v,e=e,AB=R_AB))
            if any_caps([u,v],rd):
                if rd.nodes[v]['genome']==GENOME1:
                    add_constraint(ilp,'c11',"r{Ab}_{e} + r{Aa}_{e} + y_{u} + y_{v} >= 1".format(Ab=R_Ab,e=e,Aa=R_Aa,u=u,v=v))
                    for r in [R_Ab,R_AB]:
                        add_constraint(ilp,'c12',"r{R}_{e} - b_{u} - b_{v} <= 0".format(R=r,u=u,v=v,e=e))
                    for i in [u,v]:
                        for r in [R_Ab,R_Aa]:
                            add_constraint(ilp,'c14',"y_{i} + {i} r{R}_{e} <= {i}".format(i=i,e=e,R=r))
                else:
                    add_constraint(ilp,'c11',"r{Ba}_{e} + r{Bb}_{e} + y_{u} + y_{v} >= 1".format(Ba=R_Ba, Bb=R_Bb,u=u,v=v,e=e))
                    add_constraint(ilp,'c12',"r{Ba}_{e} + b_{u} + b_{v} <= 2".format(u=u,v=v,Ba=R_Ba,e=e))
                    for i in [u,v]:
                        for r in [R_Bb,R_Ba]:
                            add_constraint(ilp,'c14',"y_{i} + {i} r{R}_{e} <= {i}".format(i=i,e=e,R=r))
            elif rd.nodes[u]['genome']==GENOME1:
                se = maybe_adjacent_selfedge(rd,u)
                se.extend(maybe_adjacent_selfedge(rd,v))
                se=list(set(se))
                if len(se) == 0:
                    add_constraint(ilp,'c13',"r{ab}_{e} = 0".format(ab=R_ab,e=e))
                else:
                    add_constraint(ilp,'c13',"{sumstr} - r{ab}_{e} >= 0".format(sumstr=' + '.join(['x_{e}'.format(e=e) for e in se]),ab=R_ab,e=e))
                

            
            


                    

                

def sibs_constraints(sibs, ilp, matching_model):
    '''
    Implement Constraints c.02, c.04 .
    '''
    for fam in sibs:
        for d_, e_ in sibs[fam]:
            d = canon_disp_tup(d_)
            e = canon_disp_tup(e_)
            #c.04
            cs = 'x_%s - x_%s = 0'%(d,e)
            add_constraint(ilp, 'c04', cs)
        #c.02
        headsum = ' + '.join(['x_%s'%canon_disp_tup(d) for d,_ in sibs[fam]])
        lower, upper = matching_model[fam]
        ch1 = '%s >= %i'%(headsum, lower)
        ch2 = '%s <= %i'%(headsum, upper)
        add_constraint(ilp,'c02',ch1)
        add_constraint(ilp,'c02',ch2)
        

MIN_CODE = 0
MAX_CODE = 8
SINGLE_LABS = ['A','B','a','b']
MINS_CODE = 0
MAXS_CODE = 16
code = dict(zip('abAB',[1,2,4,8]))

sumcode = dict([(a1+a2,c1+c2) for a1,c1 in code.items() for a2,c2 in code.items() if a1<a2]) 
sumcode['c'] = 0

def vertex_constraints(rd, ilp,bm=False):
    '''
    Implement constraints c02, c06.
    '''
    #print(sumcode)
    for v,data in rd.nodes(data=True):
        insertl(ilp, 'binary', 'z_%i'%v)
        if not bm:
            insertl(ilp,'general',('m_%i'%v,MIN_CODE,MAX_CODE))
            insertl(ilp,'general',('n_%i'%v,MIN_CODE,MAX_CODE))
        else:
            for x in SINGLE_LABS:
                insertl(ilp,'binary','m_{v}_{x}'.format(v=v,x=x))
                insertl(ilp,'binary','n_{v}_{x}'.format(v=v,x=x))
        insertl(ilp, 'general', ('y_%i'%v, 0, v))
        for xy in sumcode:
            insertl(ilp,'binary','r{xy}_{v}'.format(xy=xy,v=v))
        #c.03
        adjacent_extedges = [canon_disp_tup((v,u,k)) for u, edges in rd[v].items() for k,data in edges.items() if data['etype'] == EXTREMITYEDGE]
        adjacent_selfedges = [canon_disp_tup((v,u,k)) for u, edges in rd[v].items() for k,data in edges.items() if data['etype'] == SELFEDGE]
        smd = ' + '.join(['d_%s'%e for e in adjacent_selfedges] + ['x_%s'%e for e in adjacent_extedges])
        add_constraint(ilp, 'c03', '%s = 1'%smd)
        #c.06
        cs = '%i z_%i - y_%i <= 0'%(v,v,v)
        add_constraint(ilp, 'c06',cs)
        if data['telomere'] == True:
            #c.07
            if not bm:
                if data['genome'] == GENOME1:
                    vcode = code['A']
                elif data['genome'] == GENOME2:
                    vcode = code['B']
                cs = 'n_{v} = {vcode}'.format(v=v,vcode=vcode)
                add_constraint(ilp,'c07',cs)
            else:
                if data['genome']==GENOME1:
                    gnm = 'A'
                else:
                    gnm = 'B'
                cs = 'n_{v}_{gnm} = 1'.format(v=v,gnm=gnm)
                add_constraint(ilp,'c07',cs)
        if bm:
            #c.09
            sm = ' + '.join(['m_{v}_{x}'.format(v=v,x=x) for x in SINGLE_LABS])
            cs = '{sm} <= 1'.format(sm=sm)
            add_constraint(ilp,'c09',cs)
            sn = ' + '.join(['n_{v}_{x}'.format(v=v,x=x) for x in SINGLE_LABS])
            cs = '{sn} <= 1'.format(sn=sn)
            add_constraint(ilp,'c09',cs)
        #c.11
        if not bm:
            cs = 'm_{v} - n_{v} - 8 z_{v} <= 0'.format(v=v)
            add_constraint(ilp,'c11',cs)
            cs = 'n_{v} - m_{v} - 8 z_{v} <= 0'.format(v=v)
            add_constraint(ilp,'c11',cs)
        else:
            for lab in SINGLE_LABS:
                cs = 'm_{v}_{l} - n_{v}_{l} - z_{v} <= 0'.format(v=v,l=lab)
                add_constraint(ilp,'c11',cs)
                cs = 'n_{v}_{l} - m_{v}_{l} - z_{v} <= 0'.format(v=v,l=lab)
                add_constraint(ilp,'c11',cs)
        #c.12
        sumstr = ' + '.join(['r{xy}_{v}'.format(xy=xy,v=v) for xy in sumcode])
        cs = '{sum} - z_{v} = 0'.format(sum=sumstr,v=v)
        add_constraint(ilp,'c12',cs)
        if not bm:
            for xy, codexy in sumcode.items():
                #c.13
                cs = '{K} r{xy}_{v} - m_{v} - n_{v} <= 0'.format(v=v,xy=xy,K=codexy)
                add_constraint(ilp,'c13',cs)
                #c.14
                cs = 'm_{v} + n_{v} + 16 r{xy}_{v} <= {kplussixteen}'.format(v=v,xy=xy,kplussixteen=codexy+16)
                add_constraint(ilp,'c14',cs)
        else:
            for xy in sumcode:
                if not xy=='c':
                    x = xy[0]
                    y = xy[1]
                    for z in xy:
                        #c.13
                        cs = 'r{xy}_{v} - m_{v}_{z} - n_{v}_{z} <= 0'.format(xy=xy,v=v,z=z)
                        add_constraint(ilp,'c13',cs)
                    #c.14
                    cs = 'r{xy}_{v} - m_{v}_{x} - n_{v}_{y} >= -1'.format(v=v,x=x,y=y,xy=xy)
                    add_constraint(ilp,'c14',cs)
                    cs = 'r{xy}_{v} - m_{v}_{y} - n_{v}_{x} >= -1'.format(v=v,x=x,y=y,xy=xy)
                    add_constraint(ilp,'c14',cs)
                else:
                    #c.13
                    sm = ' + '.join(['m_{v}_{x}'.format(v=v,x=x) for x in SINGLE_LABS])
                    sn = ' + '.join(['n_{v}_{x}'.format(v=v,x=x) for x in SINGLE_LABS])
                    cs = '8 rc_{v} + {sm} + {sn} <= 8'.format(v=v,sm=sm,sn=sn)
                    add_constraint(ilp,'c13',cs)
        insertl(ilp, 'obj', (-1, 'rc_%i'%v))


def vertex_constraints_n(rd, ilp):
    '''
    Implement constraints c02, c05, c07
    '''
    #print(sumcode)
    for v,data in rd.nodes(data=True):
        insertl(ilp, 'binary', 'z_%i'%v)
        insertl(ilp,'binary','b_{v}'.format(v=v))
        insertl(ilp, 'general', ('y_%i'%v, 0, v))
        #c.02
        adjacent_extedges = [canon_disp_tup(((v,u,k))) for u, edges in rd[v].items() for k,data in edges.items() if data['etype'] == EXTREMITYEDGE]
        adjacent_selfedges = [canon_disp_tup(((v,u,k))) for u, edges in rd[v].items() for k,data in edges.items() if data['etype'] == SELFEDGE]
        if not any_caps([v],rd):
            smd = ' + '.join(['x_%s'%e for e in adjacent_selfedges] + ['x_%s'%e for e in adjacent_extedges])
            add_constraint(ilp, 'c02', '%s = 1'%smd)
        #c.05
        cs = '%i z_%i - y_%i <= 0'%(v,v,v)
        add_constraint(ilp, 'c05',cs)
        if data['telomere'] == True:
            #c.07
            if data['genome'] == GENOME1:           
                cs = 'b_{v} = 0'.format(v=v)
                add_constraint(ilp,'c07',cs)
            else:
                cs = 'b_{v} = 1'.format(v=v)
                add_constraint(ilp,'c07',cs)


def singleton_constraints(circs, ilp,dvar='d',conss='c18'):
    for i, c in enumerate(circs):
        #c.18
        k = len(c)
        smd = ' + '.join(['%s_%s'%(dvar,canon_disp_tup(e)) for e in c])
        add_constraint(ilp, conss, '%s - s_%i <= %i'%(smd, i, k-1))
        insertl(ilp, 'obj', (1, 's_%i'%i))
        insertl(ilp, 'binary', 's_%i'%i)

def summation_constraints(rd,ilp):
    #c.15,c.16
    sumstr = ' + '.join(['rAa_{v}'.format(v=v) for v in rd.nodes])
    cs = '{sum} -  rABa <= 0'.format(sum=sumstr)
    add_constraint(ilp,'c15',cs)
    sumstr = ' + '.join(['rBa_{v}'.format(v=v) for v in rd.nodes])
    cs = '{sum} -  rABa <= 0'.format(sum=sumstr)
    add_constraint(ilp,'c15',cs)
    sumstr = ' + '.join(['rAb_{v}'.format(v=v) for v in rd.nodes])
    cs = '{sum} -  rABb <= 0'.format(sum=sumstr)
    add_constraint(ilp,'c16',cs)
    sumstr = ' + '.join(['rBb_{v}'.format(v=v) for v in rd.nodes])
    cs = '{sum} -  rABb <= 0'.format(sum=sumstr)
    add_constraint(ilp,'c16',cs)
    #c.17
    transitions = ' + '.join(['rab_{v}'.format(v=v) for v in rd.nodes])
    oddpaths =    ' - '.join(['rAB_{v}'.format(v=v) for v in rd.nodes])
    cs = '{transitions} + rABa + rABb - {oddpaths} - 2 q <= 0'.format(transitions=transitions,oddpaths=oddpaths)
    add_constraint(ilp,'c17',cs)
    insertl(ilp, 'obj', (1, 'q'))
    insertl(ilp, 'general',('q',None,None))
    insertl(ilp, 'general',('rABa',0,None))
    insertl(ilp, 'general',('rABb',0,None))

def summation_constraints_n(rd,ilp):
    sumstr = ' + '.join(['0.5 x_{e}'.format(e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==EXTREMITYEDGE])
    add_constraint(ilp,'c15','{sum} - n = 0'.format(sum=sumstr))
    insertl(ilp,'general',('n',0,None))
    sumstr = ' + '.join(['r{c}_{e}'.format(c=R_C,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and not any_caps([u,v],rd)])
    add_constraint(ilp,'c16','{sum} - c = 0'.format(sum=sumstr))
    insertl(ilp,'general',('c',0,None))
    sumstr = 'piiab + ptta + pttb - pttAB - 2 q <= 0'
    add_constraint(ilp,'c17',sumstr)
    insertl(ilp, 'general',('q',None,None))
    for p in ['piiab','ptta','pttb','pttAB']:
        insertl(ilp,'general',(p,0,None))
    sumstr = ' + '.join(['r{ab}_{e}'.format(ab=R_ab,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and not any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME1])
    add_constraint(ilp,'c18','{sumstr} - piiab = 0'.format(sumstr=sumstr))
    sumstr = ' + '.join(['r{Ab}_{e}'.format(Ab=R_Ab,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME1])
    add_constraint(ilp,'c19','{sumstr} - p{Ab} = 0'.format(sumstr=sumstr,Ab=R_Ab))
    sumstr = ' + '.join(['r{Aa}_{e}'.format(Aa=R_Aa,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME1])
    add_constraint(ilp,'c21','{sumstr} - p{Aa} = 0'.format(sumstr=sumstr,Aa=R_Aa))
    sumstr = ' + '.join(['r{Ba}_{e}'.format(Ba=R_Ba,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME2])
    add_constraint(ilp,'c20','{sumstr} - p{Ba} = 0'.format(sumstr=sumstr,Ba=R_Ba))
    sumstr = ' + '.join(['r{Bb}_{e}'.format(Bb=R_Bb,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME2])
    add_constraint(ilp,'c22','{sumstr} - p{Bb} = 0'.format(sumstr=sumstr,Bb=R_Bb))
    add_constraint(ilp,'c23','ptta - p{Aa} >= 0'.format(Aa=R_Aa))
    add_constraint(ilp,'c24','ptta - p{Ba} >= 0'.format(Ba=R_Ba))
    add_constraint(ilp,'c25','pttb - p{Bb} >= 0'.format(Bb=R_Bb))
    add_constraint(ilp,'c26','pttb - p{Ab} >= 0'.format(Ab=R_Ab))
    for r in [R_Ab,R_Aa,R_Ba,R_Bb]:
        insertl(ilp,'general',('p{r}'.format(r=r),0,None))
    sumstr = ' + '.join(['r{AB}_{e}'.format(AB=R_AB,e=canon_disp_tup((u,v,k))) for u,v,k,etype in rd.edges(keys=True,data='etype') if etype==ADJACENCY and any_caps([u,v],rd) and rd.nodes[v]['genome']==GENOME1])
    add_constraint(ilp,'c27','{sumstr} - pttAB = 0'.format(sumstr=sumstr))
    insertl(ilp,'obj',(1,'n'))
    insertl(ilp,'obj',(-1,'c'))
    insertl(ilp,'obj',(1,'q'))

def component_optimization(rd,ilp,matching_model,gnms):
    amults = get_multiplicities(gnms[0])
    bmults = get_multiplicities(gnms[1])
    rd_ = deepcopy(rd)
    for u,v,k,data in rd.edges(keys=True,data=True):
        if data['etype'] == SELFEDGE:
            rd_.remove_edge(u,v,key=k)
    for c in nx.connected_components(rd_):
        has = {}
        has['a'] = False
        has['b'] = False
        has['A'] = False
        has['B'] = False
        for v in c:
            if rd_.nodes[v]['telomere']:
                if rd_.nodes[v]['genome'] == GENOME1:
                    has['A'] = True
                else:
                    has['B'] = True
            g = rd_.nodes[v]['gene']
            low_bnd = matching_model[g][0]
            if amults.get(g,0) > low_bnd:
                has['a'] = True
            if bmults.get(g,0) > low_bnd:
                has['b'] = True
        for x in has:
            if has[x]:
                continue
            for v in c:
                cs = 'm_{v}_{x} = 0'.format(v=v,x=x)
                add_constraint(ilp,'c19',cs)
                
        


def print_constraints(ilp, ctype, file=sys.stdout):
    for num, cs in getl(ilp, ctype):
        print(' %s.%i: %s'%(ctype,num, cs),file=file)

def print_all_constraints(ilp, file=sys.stdout):
    for ctype in ['c%02d'%i for i in range(1,MAX_CONSTRAINTS)]:
        print_constraints(ilp, ctype, file=file)

def print_objective(ilp, file=sys.stdout):
    obj_str = ' + '.join(['%f %s'%(k, var) for (k, var) in ilp['obj'] if k > 0])
    obj_str+=' - '
    obj_str+=' - '.join(['%f %s'%(-k, var) for (k, var) in ilp['obj'] if k < 0])
    print(' obj: %s'%obj_str, file=file)

def print_domains(ilp, file=sys.stdout):
    for var, lower, upper in ilp['general']:
        if upper is None:
            if lower is None:
                print(' %s free'%var,file=file)
            else:
                print(' %i <= %s <= Inf'%(lower,var),file=file)
        elif lower is None:
            print(' -Inf <= %s <= %i'%(var,upper))
        else:
            print(' %i <= %s <= %i'%(lower, var, upper), file=file)

def print_generals(ilp, file=sys.stdout):
    for var, _ , _ in ilp['general']:
        print(' %s'%var, file=file)

def print_binaries(ilp, file=sys.stdout):
    for var in ilp['binary']:
        print(' %s'%var, file=file)

def print_ilp(ilp, file=sys.stdout):
    print('Minimize', file=file)
    print_objective(ilp, file=file)
    print('Subject To', file=file)
    print_all_constraints(ilp, file=file)
    print('Bounds', file=file)
    print_domains(ilp, file=file)
    print('General', file=file)
    print_generals(ilp, file=file)
    print('Binary', file=file)
    print_binaries(ilp, file=file)
    print('End', file=file)
    

def main():
    parser = ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-mm','--maximal', action='store_true', help='Set matching model to maximal matching.')
    group.add_argument('-em','--exemplary', action='store_true', help='Set matching model to exemplary matching.')
    group.add_argument('-im','--intermediate', action='store_true', help='Set matching model to intermediate matching.')
    group.add_argument('-r', '--range', type=matching_range, nargs=2, help='Provide upper and lower percentiles to be matched per marker in range [0,1]. Actual discrete bounds will always be rounded up.')
    parser.add_argument('-c','--custom', type=FileType('r'), action='store', help='Provide a custom matching file.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--sparse-mode',action='store_true',
        help='Activate sparse version of the ILP. Default setting. Recommended.' if SHOW_LEGACY_VERSIONS else ap.SUPPRESS)
    group.add_argument('--slow-mode',action='store_false',dest='sparse_mode',
        help='Activate legacy ILP. Very slow. Not recommended.' if SHOW_LEGACY_VERSIONS else ap.SUPPRESS)
    group.add_argument('--binary-mode',action='store_true',dest='binary_mode',
        help='Activate legacy ILP, but in binary mode. Slow. Not recommended.' if SHOW_LEGACY_VERSIONS else ap.SUPPRESS)
    #group.add_argument('--no-binary-mode',action='store_false',dest='binary_mode')
    #parser.set_defaults(binary_mode=True)
    parser.add_argument('--optimize-comps',action='store_true',dest='optimize_comps',
        help='Activate component optimization with legacy ILP' if SHOW_LEGACY_VERSIONS else ap.SUPPRESS)
    #parser.add_argument('--no-optimize-comps',action='store_false',dest='optimize_comps')
    add_unimog_parsing_groups(parser)
    writewhat = parser.add_mutually_exclusive_group(required=True)
    writewhat.add_argument('--writemodel', type=FileType('w'), help='Write the matching model to a file in order to customize it.')
    writewhat.add_argument('--writeilp', type=FileType('w'), help='Write the resulting ILP to the specified file.')
    parser.set_defaults(sparse_mode=True)
    args = parser.parse_args()
    if args.binary_mode:
        args.sparse_mode = False
    #print(args.sparse_mode,args.binary_mode)
    if args.optimize_comps and not args.binary_mode:
        LOG.warning("--optimize-comps only works for binary mode. Ignoring...")
        args.optimize_comps=False
    if args.range:
        if min(args.range) == 0.0:
            LOG.warning('Lowest percentile to be matched is 0. If this is not counteracted by deletion penalties, almost only whole chromosome deletions will be modeled!')
    genomes = read_genomes(args)
    model = create_matching_model(args, genomes)
    if args.writemodel:
        write_matching_model(model, args.writemodel)
        sys.exit(0)
    rd, gi1, gi2, circs, exts = full_relational_diagram(genomes, LOG,capping=CAPPING_PSEUDO if args.sparse_mode else CAPPING_NONE)
    LOG.debug('Graph:')
    for c in rd.nodes(data=True):
        LOG.debug(c)
    for c in rd.edges(data=True):
        LOG.debug(c)
    sibs = siblings(rd, gi1, gi2)
    ilp = {}
    if args.sparse_mode:
        edge_constraints_n(rd, ilp)
    else:
        edge_constraints(rd,ilp,bm=args.binary_mode)
    sibs_constraints(sibs, ilp, model)
    if args.sparse_mode:
        vertex_constraints_n(rd, ilp)
    else:
        vertex_constraints(rd, ilp,bm=args.binary_mode)
    if args.sparse_mode:
        singleton_constraints(circs, ilp,dvar='x',conss='c28')
    else:
        singleton_constraints(circs, ilp)
    if args.sparse_mode:
        summation_constraints_n(rd, ilp)
    else:
        summation_constraints(rd, ilp)
    if args.optimize_comps and not args.sparse_mode:
        component_optimization(rd,ilp,model,genomes)
    print_ilp(ilp, file=args.writeilp)
    

st = logging.StreamHandler(sys.stderr)
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

st.setLevel(logging.DEBUG)
st.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
LOG.addHandler(st)


main()

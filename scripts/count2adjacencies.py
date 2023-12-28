#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from sys import stdout, stderr, exit
from collections import defaultdict
from itertools import combinations
from math import comb
from functools import reduce
from os.path import basename, dirname, join
import logging
import csv

# import from third-party packages
import pandas as pd
import numpy as np
import ete3

# import from own packages


#
# global variables
#

ORIENT_FORWARD = 1
ORIENT_REVERSE = 0

EXTR_HEAD = 'h'
EXTR_TAIL = 't'
EXTR_CAP  = 'o'

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

ids = pd.IndexSlice


def readCounts(filename):
    header_start = 0
    with open(filename) as f:
        for line in f:
            if line.startswith('#|'):
                header_start += 1
            else:
                break
    return pd.read_csv(filename, sep='\t', index_col=0, header=[header_start], low_memory=False)


def readGF(data):
    _df = pd.read_csv(data, sep='\t', header=None, index_col=0, names=['Members'])
    _df['count'] = _df.Members.map(lambda x: len(x.split()))
    df = pd.DataFrame(data = {'family': reduce(lambda a, b: a + b, map(lambda x: [x[0]] * x[1]['count'], _df.iterrows())),
                                 'genome_gene': reduce(lambda a, b: a+b, _df.Members.map(lambda x: x.split()))
                                 })
    df['genome'] = df.genome_gene.map(lambda x: x.split('|', 1)[0])
    # intentially programmed to fail in case | is not part of string...
    df['gene'] = df.genome_gene.map(lambda x: x.split('|', 1)[1])
    del df['genome_gene']
    df.sort_values(['genome', 'gene'], inplace=True)
    df.set_index(['genome', 'gene'], inplace=True)

    return df


def readGO(data, path_prefix):
    df = pd.DataFrame({
        'genome_gene': pd.Series(dtype=str),
        'orientation': pd.Series(dtype=int),
        'start': pd.Series(dtype=int),
        'end': pd.Series(dtype=int),
        'chromosome': pd.Series(dtype=str),
        })
    for genome, f in csv.reader(data, delimiter='\t'):
        fname = join(path_prefix, f)
        LOG.info('++ parsing gene order of {} from {}'.format(genome, fname))
        df = pd.concat([df,
                           pd.read_csv(fname, sep='\t', header=None,
                                       names=['genome_gene', 'orientation', 'start', 'end', '_unused_', 'chromosome'])],
                          join='inner', ignore_index=True)

    df['genome'] = df.genome_gene.map(lambda x: x.split('|', 1)[0])
    # intentially programmed to fail in case | is not part of string...
    df['gene'] = df.genome_gene.map(lambda x: x.split('|', 1)[1])
    df.sort_values(['genome', 'chromosome', 'start', 'end'], inplace=True)
    df.set_index(['genome', 'gene'], inplace=True)

    # remove columns that we no longer need
    del df['genome_gene']
    del df['start']
    del df['end']

    return df


def constructExtantAdjacenciesTable(df_go):
    """ Constructs adjacency table from given gene order table. """

    map_orient1 = {0: EXTR_TAIL, 1: EXTR_HEAD, 2: EXTR_CAP}
    map_orient2 = {0: EXTR_HEAD, 1: EXTR_TAIL, 2: EXTR_CAP}

    df = pd.DataFrame({
        'family1': pd.Series(dtype=str),
        'gene1': pd.Series(dtype=str),
        'ext1': pd.Series(dtype=str),
        'family2': pd.Series(dtype=str),
        'gene2': pd.Series(dtype=str),
        'ext2': pd.Series(dtype=str),
        })
    # at this point, we assume that each extant chromosome is linear
    for chrom in df_go.chromosome.unique():
        # the following code creates an "adjacency table" from df_go joining it with a copy of itself that is shifted by one; then telomeres are
        # added and the table is prepared so that columns match that of the resulting df table
        df_c = df_go.loc[df_go.chromosome == chrom]
        df_c1 = df_c.shift(1) 
        df_c1.loc[0, ['gene', 'orientation', 'family']] = ['0', 2, 't']
        df_c12 = df_c1.join(df_c, lsuffix='1', rsuffix='2')
        df_c12['ext1'] = df_c12.orientation1.map(map_orient1.get)
        df_c12['ext2'] = df_c12.orientation2.map(map_orient2.get)
        last = df_c12.tail(1)
        cols = ['gene1', 'family1', 'ext1', 'gene2', 'family2', 'ext2']
        df_c12.loc[last.index.item() + 1, cols] = [last.gene2.item(), last.family2.item(), map_orient1[last.orientation2.item()],
                                                   str(last.index.item()+1), 't', EXTR_CAP]
        canonizeAdjacencies(df_c12)
        df = pd.concat([df , df_c12], join='inner', ignore_index=True)

    return df


def canonizeAdjacencies(df):
    """ makes sure adjacencies are represented in canonical form """
    sel_unsrtd = df.apply(lambda x: (x.family1, x.ext1, x.gene1) > (x.family2, x.ext2, x.gene2), axis=1)
    df_tmp1 = df.loc[sel_unsrtd, ['family1', 'ext1', 'gene1']]
    df_tmp1.columns = ['family2', 'ext2', 'gene2']
    df_tmp2 = df.loc[sel_unsrtd, ['family2', 'ext2', 'gene2']]
    df_tmp2.columns = ['family1', 'ext1', 'gene1']
    df.loc[sel_unsrtd, ['family1', 'ext1', 'gene1']] = df_tmp2
    df.loc[sel_unsrtd, ['family2', 'ext2', 'gene2']] = df_tmp1


def constructExtantAdjacenciesTableAll(tree, df_go):
    tc = 0
    df = pd.DataFrame({
        'species': pd.Series(dtype=str),
        'family1': pd.Series(dtype=str),
        'gene1': pd.Series(dtype=str),
        'ext1': pd.Series(dtype=str),
        'family2': pd.Series(dtype=str),
        'gene2': pd.Series(dtype=str),
        'ext2': pd.Series(dtype=str),
        })
    # construct table of all observed adjacencies in extant genomes
    for v in tree.get_leaves():
        df_v = constructExtantAdjacenciesTable(df_go.loc[ids[v.name,:,:]].reset_index())
        df_v['species'] = v.name
        df = pd.concat([df, df_v], ignore_index=True)
    # initialize weight of extant adjacencies with 1
    df['weight'] = 1.0
    return df


def recruitAncestralAdjacencies(tree, df_extant):
    """
    Assumes that genes in df_go are already ordered by their genomic position.

    This function performs a bottom-up/top-down swipe to construct a weighted adjacency set for each internal node of the given phylogeny. The result
    will be such that all internal nodes will have the same adjacencies, but their counts will be different, depending on how many times the
    corresponding adjacency has been seen on any path through the corresponding node.

    Note that the code may sometimes look weird, because it works on general trees, but typically we expect the tree to be binary.
    """
    # ignore gene associations of extant genes for bottom-up traversal and ignore duplicate adjacencies
    df = df_extant[['species', 'family1', 'ext1', 'family2', 'ext2', 'weight']].drop_duplicates()

    #
    # initialize table for all ancestral genomes 
    #
    df_template = df[['family1', 'ext1', 'family2', 'ext2']].drop_duplicates().set_index(['family1', 'ext1', 'family2', 'ext2'])
    df_template['weight'] = 0
    df.set_index(['species', 'family1', 'ext1', 'family2', 'ext2'], inplace=True)
    for v in tree.traverse():
        df = pd.concat([df, pd.concat({v.name: df_template}, names=['species'])])
    # we created some duplicates for leaves, so we have to drop those
    df = df[~df.index.duplicated(keep='first')]
    df.sort_index(inplace=True)

    return df.loc[ids[[x.name for x in tree.traverse() if not x.is_leaf()], :, :, :]].reset_index()


def weighAdjacenciesByPathConservation(df_anc, df_extant, tree):

    # introduce temporary columns "up_count" and "down_count"
    df = pd.concat([df_extant[['species', 'family1', 'ext1', 'family2', 'ext2', 'weight']].drop_duplicates(), df_anc], ignore_index=True)
    df['up_count'] = df['weight']
    df['down_count'] = 0
    df.set_index(['species', 'family1', 'ext1', 'family2', 'ext2'], inplace=True)
    df.sort_index(inplace=True)

    #
    # bottom-up traversal
    #
    for v in tree.traverse('postorder'):
        for u, w in combinations(v.get_children(), 2):
            df_u = pd.concat({v.name: df.loc[ids[u.name,:, :, :, :]]}, names=['species'])
            df_w = pd.concat({v.name: df.loc[ids[w.name,:, :, :, :]]}, names=['species'])
            # set weight, which is the product of the children's counts
            # by construction, df_u and df_w have same index
            df.loc[ids[v.name, :, :, :], 'weight'] += df_u.up_count * df_w.up_count

        # update count of v to the sum of the children's counts
        for u in v.get_children():
            df_u = pd.concat({v.name: df.loc[ids[u.name,:, :, :, :]]}, names=['species'])
            df.loc[df_u.index, 'up_count'] += df_u.up_count
    #
    # top-down traversal
    # 
    for v in tree.traverse('preorder'):
        df.loc[ids[v.name, :, :, :], 'weight'] += df.loc[ids[v.name, :, :, :], 'up_count'] * df.loc[ids[v.name, :, :, :], 'down_count']
        # update count for siblings
        df_children = df.loc[ids[[u.name for u in v.get_children()], :, :, :]].groupby(['family1', 'ext1', 'family2', 'ext2']).sum()
        for u in v.get_children():
            s = pd.concat({u.name: df_children.down_count + df.loc[v.name, 'up_count']}, names=['species'])
            df.loc[ids[u.name, :, :, :], 'down_count'] = s - df.loc[ids[u.name, :, :, :], 'up_count']

    # we report *only* the weights of ancestral adjacencies
    df_paths = countPaths(speciesTree)
    return normalizeWeights(df_paths, df.loc[ids[[x.name for x in tree.traverse() if not x.is_leaf()], :, :, :]].reset_index())



def countPaths(tree):
    """ reports the number of leaf-to-leaf path that go through each node of the tree """

    df = pd.DataFrame(index=map(lambda x: x.name, tree.traverse()), data=0, columns=['paths', 'up_leaves', 'down_leaves'])

    # bottom-up traversal
    for v in tree.traverse('postorder'):
        if v.is_leaf():
            df.loc[v.name, 'up_leaves'] = 1
        else:
            for u, w in combinations(v.get_children(), 2):
                df.loc[v.name, 'paths'] += df.loc[[u.name, w.name], 'up_leaves'].product()
            # update leaf counter
            df.loc[v.name, 'up_leaves'] = df.loc[[u.name for u in v.get_children()], 'up_leaves'].sum()

    # top-down traversal
    for v in tree.traverse('preorder'):
        df.loc[v.name, 'paths'] += df.loc[v.name, 'up_leaves'] * df.loc[v.name, 'down_leaves']

        # update leaf count for children
        c = df.loc[[u.name for u in v.get_children()], 'up_leaves'].sum()
        for u in v.get_children():
            df.loc[u.name, 'down_leaves'] = df.loc[v.name, 'down_leaves'] + c - df.loc[u.name, 'up_leaves']

    return df[['paths']]

def normalizeWeights(df_paths, df_adjs):
    """ normalizes ancestral adjacencies """

    for s in df_adjs['species'].unique():
        df_adjs.loc[df_adjs.species == s, 'weight'] /= df_paths.loc[s, 'paths']

    return df_adjs

def weighAdjacenciesByWeightScheme(df_extant, df_anc, df_ws):

    ids = pd.IndexSlice
    adjs = set(df_anc[['family1', 'ext1', 'family2', 'ext2']].itertuples(index=False))
    df_anc.set_index(['species', 'family1', 'ext1', 'family2', 'ext2'], inplace=True)
    df_extant = df_extant.set_index(['family1', 'ext1', 'family2', 'ext2'])
    df_extant.sort_index(inplace=True)

    for family1, ext1, family2, ext2 in adjs:
        sp = set(df_extant.loc[(family1, ext1, family2, ext2), 'species'])
        df_w = df_ws.loc[tuple((x[1] in sp and 1 or 0 for x in df_ws.index.names))]
        df_w.index = df_w.index.droplevel(0)
        df_anc.loc[ids[df_w.index, family1, ext1, family2, ext2], 'weight'] = df_w.array

    df_anc.reset_index(inplace=True)
    return df_anc


def instantiateGenes(df_adjs, df_counts):

    df = pd.DataFrame({
        'species': pd.Series(dtype=str),
        'family1': pd.Series(dtype=str),
        'gene1': pd.Series(dtype=str),
        'ext1': pd.Series(dtype=str),
        'family2': pd.Series(dtype=str),
        'gene2': pd.Series(dtype=str),
        'ext2': pd.Series(dtype=str),
        'weight': pd.Series(dtype=int)
        })

    for anc in df_counts.columns:
        df_genes = pd.DataFrame(data = {'family': reduce(lambda a,b: a+b, map(lambda x: [x[0]] * x[1], df_counts[anc].items()))})
        df_genes.index.name = 'gene'
        df_genes.reset_index(inplace=True)
        df_genes.gene = df_genes.gene.astype(str)

        df_anc = df_adjs.loc[df_adjs.species==anc].set_index('family1').join(df_genes.set_index('family'), how='inner').reset_index()
        df_anc = df_anc.set_index('family2').join(df_genes.set_index('family'), how='inner', lsuffix='1', rsuffix='2').reset_index()
        canonizeAdjacencies(df_anc)
        df_anc.drop_duplicates(inplace=True)
        df = pd.concat([df, df_anc], ignore_index=True)

    return df


def cleanupAncestralAdjacencies(df_adjs, df_counts):

    # remove adjacencies between same extremities
    sel_rm_same_ext = (df_adjs.gene1 != df_adjs.gene2) | (df_adjs.ext1 != df_adjs.ext2)

    # prohibit singleton genes by removing selfloops from duplicates
    sel_is_dup = df_adjs.apply(lambda x: df_counts.loc[x.family1, x.species] > 1, axis=1)
    sel_rm_selfloop = ~sel_is_dup | (df_adjs.gene1 != df_adjs.gene2)

    return df_adjs.loc[sel_rm_same_ext & sel_rm_selfloop]


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=str,
            help='phylogenetic tree in newick format')
    parser.add_argument('gf_counts', type=open,
            help='gene family count table from Miklós Csűrös\' Count software')
    parser.add_argument('gene_families', type=open,
            help='gene-to-family assignment table')
    parser.add_argument('gene_orders', type=open,
            help='file pointing to gene order tables of extant species')
    parser.add_argument('-p', '--use_weight_scheme', type=open,
            help='use weight specified as a function of leaf absence/presence in given file')
    parser.add_argument('-m', '--min_weight', type=float, default=0.0,
            help='minimum weight of ancestral adjacencies')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    # load & process input data
    LOG.info('loading species tree from {}'.format(args.tree))
    speciesTree = ete3.Tree(args.tree, format=1)

    LOG.info('loading gene family counts from {}'.format(args.gf_counts.name))
    df_counts = readCounts(args.gf_counts.name)
    LOG.info('loading gene family assignments from {}'.format(args.gene_families.name))
    df_gf = readGF(args.gene_families)

    LOG.info('loading extant gene orders from {}'.format(args.gene_orders.name))
    # read gene order table and join with df_gf table
    df_go = readGO(args.gene_orders, dirname(args.gene_orders.name)).join(df_gf)

    # let's do it!

    df_extant = constructExtantAdjacenciesTableAll(speciesTree, df_go)
    df_anc = recruitAncestralAdjacencies(speciesTree, df_extant)

    if args.use_weight_scheme:
        df_ws = pd.read_csv(args.use_weight_scheme, sep='\t', header=[0, 1])
        df_ws.set_index([x for x in df_ws.columns if x[0] == 'configuration'], inplace=True)
        df_anc = weighAdjacenciesByWeightScheme(df_extant, df_anc, df_ws)
    else:
        df_anc = weighAdjacenciesByPathConservation(df_anc, df_extant, speciesTree)

    # join extant and ancestral adjacency sets
    df = pd.concat([cleanupAncestralAdjacencies(instantiateGenes(df_anc, df_counts[df_anc.species.unique()]), df_counts), df_extant],
                   ignore_index=True)

    df.gene1 = df.apply(lambda x: '_'.join((x.family1, x.gene1)), axis=1)
    df.gene2 = df.apply(lambda x: '_'.join((x.family2, x.gene2)), axis=1)

    # enforce minimum weight
    df = df.loc[df.weight >= args.min_weight]

    # output final adjacency set
    df.rename(columns={'species': '#species'}, inplace=True)
    df[['#species', 'gene1', 'ext1', 'gene2', 'ext2', 'weight']].to_csv(stdout, sep='\t', index=False)

    LOG.info('DONE')


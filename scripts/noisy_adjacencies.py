#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
from collections import defaultdict
import random
import csv

# import from own packages
import data_utils as du

NOISY_ADJ_WEIGHT = 1
ADVERSARIAL_ADJ_WEIGHT = 1

MAX_TRIES = 10000 # Max number of tries to find a random adjacency

def createRandomAdjacency(species,genesList,adjacenciesList):
    repeat = 0
    while repeat <= MAX_TRIES:
        [gene1,gene2] = random.sample(genesList[species],2)
        sign1 = random.choice(['h','t'])
        sign2 = random.choice(['h','t'])
        ext1 = (gene1,sign1)
        ext2 = (gene2,sign2)
        if ext1>ext2:
            ext1,ext2=ext2,ext1
        if [ext1,ext2] not in adjacenciesList[species]:
            return([ext1,ext2])
        repeat+=1

def createRandomAdversarialAdjacency(species,parent,familiesDict,adjacenciesList,noiseFamilies):
    repeat = 0
    familiesAvailable = list(familiesDict[species].keys())
    while repeat <= MAX_TRIES:
        [fam1,fam2] = random.choice(noiseFamilies[parent])
        if fam1 in familiesAvailable and fam2 in familiesAvailable:
            gene1 = random.choice(familiesDict[species][fam1])
            gene2 = random.choice(familiesDict[species][fam2])
            sign1 = random.choice(['h','t'])
            sign2 = random.choice(['h','t'])
            ext1 = (gene1,sign1)
            ext2 = (gene2,sign2)
            if ext1 != ext2:
                if ext1>ext2:
                    ext1,ext2=ext2,ext1
                if [ext1,ext2] not in adjacenciesList[species]:
                    return([ext1,ext2])
            repeat+=1


def addNoiseFamily(noiseFamilies,ext1,ext2):
    fam1 = du.getFamily(ext1)
    fam2 = du.getFamily(ext2)
    if not [fam1,fam2] in noiseFamilies:
        noiseFamilies.append([fam1,fam2])

def DFS(childrenDict,root,nodesOrder):
    nodesOrder.append(root)
    for child in childrenDict[root]:
        DFS(childrenDict,child,nodesOrder)

def addAdversarialNoise(dataSet, speciesTree, nbNovelAdjacencies, probaOfRepeat, noiseInLeaves=False,
        leavesDict={}):
    '''Add noise, i.e. novel adjacencies, to a dataset Adds to each species
    nbNovelAdjacencies The adjacencies are added at random with probability 1-probaOfRepeat
    but at the root (proba 1) otherwise they are added between two genes whose families
    already form together a noisy adjacency at the previous node

    @returns a novel dataset composed of a species list, a genes list and an
    adjacency list'''
    speciesList     = dataSet['species'].copy()
    genesList       = dataSet['genes'].copy()
    adjacenciesList = dataSet['adjacencies'].copy()
    weightsDict     = dataSet['weights'].copy()
    familiesDict    = dataSet['families'].copy()

    parentDict   = {}
    childrenDict = defaultdict(list)
    for [sp1,sp2] in speciesTree:
        parentDict[sp1] = sp2
        childrenDict[sp2].append(sp1)
    speciesDFSList = []
    DFS(childrenDict,'Root',speciesDFSList)

    noiseFamilies = {species:[] for species in speciesList}

    # Random adjacencies at the root
    for i in range(nbNovelAdjacencies):
        [ext1,ext2] = createRandomAdjacency('Root',genesList,adjacenciesList)
        adjacenciesList['Root'].append([ext1,ext2])
        weightsDict['Root'][ext1,ext2] = NOISY_ADJ_WEIGHT
        addNoiseFamily(noiseFamilies['Root'],ext1,ext2)

    # Other species
    for species in speciesDFSList:
        if (not species=='Root') and \
           ((not leavesDict[species]) or (noiseInLeaves and leavesDict[species])):
            parent = parentDict[species]
            for i in range(nbNovelAdjacencies):
                adjType = random.random()
                if adjType < probaOfRepeat:
                    adj = createRandomAdversarialAdjacency(species,parent,familiesDict,adjacenciesList,noiseFamilies)
                    if adj != None:
                        ext1, ext2 = adj
                        weightsDict[species][ext1, ext2] = ADVERSARIAL_ADJ_WEIGHT
                    else:
                        print('WARNING: no suitable adjacency for ' + \
                              'adversarial noise for species '+species +\
                              ' with parent '+parent+ \
                              ' having '+str(len(noiseFamilies[parent])) +\
                              ' available noisy adjacencies', file = stderr)
                        [ext1,ext2] = createRandomAdjacency(species,genesList,adjacenciesList)
                        weightsDict[species][ext1,ext2] = NOISY_ADJ_WEIGHT
                else:
                    [ext1,ext2] = createRandomAdjacency(species,genesList,adjacenciesList)
                    weightsDict[species][ext1,ext2] = NOISY_ADJ_WEIGHT
                adjacenciesList[species].append([ext1,ext2])
                addNoiseFamily(noiseFamilies[species],ext1,ext2)

    return {'species':speciesList, 'genes':genesList,
            'adjacencies':adjacenciesList, 'weights':weightsDict,
            'families':familiesDict}


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=open,
                        help='phylogenetic tree as parent-child relation table')
    parser.add_argument('trueAdjacencies', type=open,
                        help='true adjacencies of the genomes in the phylogeny')
    parser.add_argument('noiseLevel', type=int,
                        help='number of noisy adjacencies to add to the true gene order')
    parser.add_argument('adversarialNoiseLevel', type=float,
                        help='ratio of adversarial noisy adjacencies to add to the true gene order')
    parser.add_argument('-l', '--noiseInLeaves', action = 'store_true',
                        help='adds noise to the leaves')
    parser.add_argument('outputName', type=FileType('w'),
                        help='name for the output adjacencies file')

    args = parser.parse_args()

    # load data
    speciesTree     = du.parseTree(args.tree)
    leavesDict      = du.getLeaves(speciesTree)
    trueAdjacencies = du.parseAdjacencies(args.trueAdjacencies)

    # construct noisy adjacency set
    noisyDataSet = addAdversarialNoise(trueAdjacencies, speciesTree, args.noiseLevel,
                                       args.adversarialNoiseLevel, noiseInLeaves=args.noiseInLeaves,
                                       leavesDict=leavesDict)

    # write output
    du.writeAdjacencies(noisyDataSet['adjacencies'], noisyDataSet['weights'],
            args.outputName)

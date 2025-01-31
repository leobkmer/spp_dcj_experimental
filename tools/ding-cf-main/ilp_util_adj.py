CHR_CIRCULAR = ')'
CHR_LINEAR = '|'
ORIENT_POSITIVE = '+'
ORIENT_NEGATIVE = '-'
EXTREMITY_TAIL = 't'
EXTREMITY_HEAD = 'h'
EXTREMITY_TELOMERE = '$'
#all telomeres are equivalent concerning the distance
TELOMERE_ID = 'TL'

#TODO: protect for parallel execution
class Simple_Id_Generator:
    last = 0
    def get_new(self):
        self.last+=1
        return self.last

#util for dict with lists
def insertl(d, k, v):
    if k not in d:
        d[k] = []
    d[k].append(v)

def getl(d,k):
    if not k in d:
        return []
    return d[k]

def dadd(d,k,v):
    if k not in d:
        d[k] = 0
    d[k]+=v

def dgetc(d,k):
    if not k in d:
        return 0
    return d[k]

def readGenomes(data, genomesOnly=None):
    """Read genome in UniMoG format
    (https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual)"""

    res = list()

    # helper function for parsing each individual gene
    str2gene = lambda x: x.startswith(ORIENT_NEGATIVE) and (ORIENT_NEGATIVE, \
            x[1:]) or (ORIENT_POSITIVE, x.lstrip(ORIENT_POSITIVE))
    # process each line, assuming that the file is well-formatted
    skip = False
    for line in data:
        line = line.strip()
        if line:
            if line.startswith('>'):
                genomeName = line[1:].strip()
                if genomesOnly == None or genomeName in genomesOnly:
                    skip = False
                    res.append((genomeName, list()))
                elif genomesOnly:
                    skip = True
            elif line[-1] not in (CHR_CIRCULAR, CHR_LINEAR):
                raise Exception('Invalid format, expected chromosome to ' + \
                        'end with either \'%s\' or \'%s\'' %(CHR_CIRCULAR, \
                        CHR_LINEAR))
            elif not skip:
                res[-1][1].append((line[-1], list(map(str2gene, line[:-1].split()))))
    return res

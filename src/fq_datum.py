
import pyfastx

import numpy as np
from collections import Counter

def readGCcontent(seq):
    seq = seq.upper()
    A = seq.count("A")
    T = seq.count("T")
    C = seq.count("C")
    G = seq.count("G")
    B = [A, T, C, G]
    GC = round( (C+G)/sum(B), 3)
    return GC

def readAvgQscore(quali):
    read_length = len(quali)
    value_list = Counter( list(quali) )
    value_sum = sum([ k * v for k, v in value_list.items() ])
    read_qscore = value_sum / read_length
    return read_qscore

def endsBaseDataCapture(seq, qual, shift_length):
    a = 0

def totBaseInfoCapture(seq, qual):
    b = 0

def random_readnum(seed_num, read_size, sample_num):
    np.random.seed(seed_num)
    sample_list = np.random.randint(1, read_size, (1, sample_num) )
    return sample_list

fq = pyfastx.Fastq('../test/ecoli.fq.gz')

def readParse(read, seqdict):
    seqdict['ID'].append(read.name )
    seqdict['GC'].append( readGCcontent( read.seq ) )
    seqdict['LEN'].append( len( read.seq ) )
    seqdict['QUAL1'].append( float(read.name.split('_')[-1]) )
    seqdict['QUAL2'].append( readAvgQscore( read.quali ) )
    return seqdict

def sampling_analyser(fq, seed_num, read_num):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2':[]} ) # QUAL1: read basecall Q, QUAL2: read average Q
    sample_list = random_readnum(seed_num, len(fq), read_num)
    for i in sample_list[0]:
        read = fq[i]
        seqdict = readParse(read, seqdict)
    return seqdict
def overall_analyser(fq):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2': []} )  # QUAL1: read basecall Q, QUAL2: read average Q
    for read in fq:
        seqdict = readParse(read, seqdict)
    return seqdict
seqdict = overall_analyser(fq)
print(seqdict)
import numpy as np
import pandas as pd
import math
import umap
import hdbscan
import sklearn.cluster as cluster
from numba import jit
import re
import networkx as nx
import pygraphviz
import random
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt

import scipy.sparse as sp
def document_frequency(X):
    """Count the number of non-zero values for each feature in sparse X."""
    if sp.isspmatrix_csr(X):
        return np.bincount(X.indices, minlength=X.shape[1])
    else:
        return np.diff(X.indptr)
    


def display_graph(G):
    plt.figure(1, figsize=(8, 8))
    # layout graphs with positions using graphviz neato
    pos = graphviz_layout(G, prog="neato")
    # color nodes the same in each connected subgraph
    C = (G.subgraph(c) for c in nx.connected_components(G.to_undirected()))
    for g in C:
        c = [random.random()] * nx.number_of_nodes(g)  # random color...
        nx.draw(g,
                pos,
                node_size=20,
                node_color=c,
                vmin=0.0,
                vmax=1.0,
                with_labels=False
               )
    plt.show()


def assembly(graph):
    R = nx.dag_longest_path(nx.DiGraph(graph), weight='weight', default_weight=1)
    return R[0] + ''.join([x[-1] for x in R[1:]])

def assembly_complete(graph):
    S = [nx.DiGraph(graph.subgraph(c).copy()) for c in nx.connected_components(graph.to_undirected())]
    return np.array([assembly(subg) for subg in S])

import networkx as nx



#### get kmers edge list 
def tuple_to_graph_edge(graph, kmer1, kmer2):
    if graph.has_edge(kmer1, kmer2):  graph[kmer1][kmer2]['weight'] += 1
    else: graph.add_edge(kmer1, kmer2, weight=1)
    return

### add edges from a read in form of string
def add_read_to_graph(graph, read, kmer_len):
    for i in range(1, len(read)-kmer_len+1):
        tuple_to_graph_edge(graph,read[i-1:i-1+kmer_len], read[i:i+kmer_len]) 
    return 

### add edges from a read in form of integer
def add_read_to_graph_as_int(graph, read, kmer_len):
    for i in range(1, len(read)-kmer_len+1):
        tuple_to_graph_edge(graph,convertdatatype_undirected(read[i-1:i-1+kmer_len]), convertdatatype_undirected(read[i:i+kmer_len])) 
    return

def read2kmers_int_unidirect(read, kmer_len):
    return list(map(convertdatatype_undirected, get_kmers(read, kmer_len)))

def read2EncodedRead(read, minik = 4):
    return ''.join(map(nucKmer2encodedkmer, re.findall(r'[ATCNG]{'+str(minik)+'}',read )))

def DecodeChrRead2Read(encread):
    return ''.join([encodedKmer2nuc(x) for x in encread])


def convertdatatype_undirected(kmer):
    return int(kmer.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4').replace('N', '0'), 5)


def nucKmer2encodedkmer(kmer):
    return chr(int(kmer.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4').replace('N', '0'), 5))

    
def encodedKmer2nuc(char):    
    return np.base_repr(ord(char),5).replace('1', 'A').replace('2', 'C').replace('3', 'G').replace('4', 'T').replace('0', 'N')


def int2nuc(number):
    return np.base_repr(number, base=5).replace('1', 'A').replace('2', 'C').replace('3', 'G').replace('4', 'T').replace('0', 'N')

def nuc_to_num(read):
    return read.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4').replace('N', '0')


def nuc_to_base32(read):
    return np.base_repr(int(nuc_to_num(read), 5), base=32)

def nuc_to_base32(read):
    return np.base_repr(int(nuc_to_num(read), 5), base=32)

def base32_to_nuc(val):
    return np.base_repr(int(val, 32), base=5).replace('1', 'A').replace('2', 'C').replace('3', 'G').replace('4', 'T').replace('0', 'N')



def convertdatatype(read):
    try: inverse = int(nuc_to_num(read[::-1]), 5)
    except: return -1
    verse = int(nuc_to_num(read), 5)
    if verse < inverse: return verse
    return inverse


def convertdatatype_undirected(kmer):
    try: return int(nuc_to_num(read), 5)
    except: return -1




def get_kmers(read, ksize):
    return [read[i:i+ksize] for i in range(0, len(read)-ksize+1)]


def get_kmers_from_read(read, kmer_len):
    kmers = np.array(list(map(nuc_to_base32, get_kmers(read, kmer_len))))
    return kmers[kmers != '-1']


def get_kmers_from_read_int(read, kmer_len):
    kmers = np.array(list(map(convertdatatype, get_kmers(read, kmer_len))))
    return kmers[kmers != '-1']


 
def Fastq2KmersNPY(filename, k = 30):
    with open(filename) as FQfile:
        with open(filename.split('.')[0] + '_kmers.npy', 'wb') as to_fil:
            while 1:
                if not FQfile.readline(): break
                np.save(to_fil, get_kmers_from_read(FQfile.readline(),k), allow_pickle=True)
                FQfile.readline()
                FQfile.readline()                
                

                
def fastq2kmerstxt(filename, k = 30):
    with open(filename) as FQfile:
        with open(filename.split('.')[0] + '_kmers.txt', 'w') as to_fil:
            while 1:
                if not FQfile.readline(): break
                to_fil.write(' '.join(get_kmers(FQfile.readline(), k)))
                FQfile.readline()
                FQfile.readline()


def Read_NPY_file(filename):
    stacks = np.array([])
    with open(filename,'rb') as f:
        for i in range(int(int(subprocess.check_output('wc -l ' + filename.split('_kmers')[0] + '.fastq', shell=True).split()[0])/4)):
            stacks = np.hstack([np.load(f, allow_pickle=True), stacks])
    return stacks

def fastq2kmersH5(filename, ksize = 30, enc = 'integer'):
    df = pd.HDFStore(filename.split('.')[0] + "_kmers.h5", 'w')
    if enc == 'integer':
        with open(filename) as FQfile:
            while 1:
                line = FQfile.readline()
                if not line: break
                df.append('seq',pd.DataFrame(get_kmers_from_read(FQfile.readline(), ksize), columns = ['kmers']), index=False )
                FQfile.readline()
                FQfile.readline()
    elif enc == 'string':
        with open(filename) as FQfile:
            while 1:
                line = FQfile.readline()
                if not line: break
                df.append('seq',pd.DataFrame(get_kmers(FQfile.readline(), ksize), columns = ['kmers']), index=False )
                FQfile.readline()
                FQfile.readline()
    df.close()
    
    
def fastq2tsv(filename):
    with open(filename) as FQfile:
        with open(filename.split('.')[0] + '.tsv', 'w') as to_fil:
            while 1:
                if not FQfile.readline(): break
                to_fil.write('\t'.join([FQfile.readline().replace('\n', ''), FQfile.readline().split(' ')[0][1:], FQfile.readline().replace('\n', '')]))
                to_fil.write('\n')


                
                
###### set kmer length
##### analyzers
import numba
from numba import njit


                
                
@njit 
def numbard(dataset, k):
    window = np.zeros(k)
    ret = []
    is_line = 1
    cnt = 0
    for char in dataset:
        if char == 10: 
            is_line += 1
            if is_line > 4: is_line= 1
            cnt = 0
        elif is_line == 2: 
            if char == 78: cnt = 0
            else:
                window[cnt] = char
                cnt += 1
                if cnt == k: 
                    ret.append(window)
                    cnt -= 1
                    window[:19] = window [1:20]
    return ret
        
def ngram_numba(filename, klen = 20):
    return numbard(np.fromfile(filename, dtype='int8'), klen)




###### separating for grid search
#grid = {'vect__ngram_range': [(1, 1)],
#        'tfidf__norm': ['l1', 'l2'],
#        'clf__alpha': [1e-3, 1e-4, 1e-5]}
#partial_fit_classifiers = {
#    'SGD': SGDClassifier(max_iter=5),
#    'Perceptron': Perceptron(),
#    'NB Multinomial': MultinomialNB(alpha=0.01),
#    'Passive-Aggressive': PassiveAggressiveClassifier(),
#}
from functions import *
import re
import networkx as nx
from dask.distributed import Client
from nltk.util import ngrams
from dask import delayed
from dask import bag as db
from dask import compute
from nltk.util import ngrams
from sklearn.feature_extraction.text import HashingVectorizer , TfidfVectorizer, CountVectorizer, TfidfTransformer
from sklearn.linear_model import SGDClassifier
from sklearn.pipeline import Pipeline
import dask_ml.feature_extraction.text as dask_text
import dask
from itertools import chain
from dask import bag as db
from dask import dataframe as dd
from itertools import chain
import pandas as pd
import dask.array as da
import itertools
import dask
from dask import bag as db
from itertools import chain
Kmer_length_func = 20




cdef ng(str read, int n):
    return ngrams(read, n)

tab =  str.maketrans("ACTG", "TGAC")
cdef rev_comp(seq):
    return seq.translate(tab)[::-1]

cdef char_transform(char char):
    if char == b'A': return 'T'
    if char == b'C': return 'G'
    if char == b'G': return 'C'
    if char == b'T': return 'A'
    return 'N'


cdef char_transform_py(str char):
    if char == 'A': return 'T'
    if char == 'C': return 'G'
    if char == 'G': return 'C'
    if char == 'T': return 'A'
    return 'N'

def ord_kmer_v2(str kmer, int val):
    cdef str to_ret = ''
    cdef int i = 0
    cdef str inv = ''
    while i < val:
        inv = char_transform(kmer[val - i-1])
        if kmer[i]< inv:
            return kmer
        elif kmer[i] > inv:
            return rev_comp(kmer)
        i +=1
    return kmer


def nltk_ngram_yield(text, n):
    for i in ng(text, n):
        yield ord_kmer_tup(i, n)


def tokenize_array( str text, int klen = 20):
    cdef str window  =  ''
    cdef list returnable = []
    cdef int counter = 0
    cdef int check_if_read = 1
    for letter in text:
        if  (letter !=  'A') & (letter != 'C') & (letter != 'G') & (letter != 'T'): 
            if letter == '\n': 
                check_if_read = check_if_read%4 + 1
            window = ''
            counter = 0
        elif check_if_read == 2:
            window = window + letter
            counter += 1
            if counter == klen:
                returnable =  returnable + [window] #
                window = window[1:]
                counter -= 1
    return returnable

tab =  str.maketrans("ACTG", "TGAC")
cdef rev_comp(seq):
    return seq.translate(tab)[::-1]

cdef ord_kmer(str kmer):
    cdef str kmer2 = rev_comp(kmer)
    for i, j in zip(kmer, kmer2):
        if ord(i)< ord(j):
            return kmer
        elif ord(j)< ord(i):
            return kmer2
    return kmer

def tokenize_read( str text, int klen = 30):
    cdef str window  =  ''
    cdef list returnable = []
    cdef int counter =0
    for letter in text:
        if  (letter !=  'A') & (letter != 'C') & (letter != 'G') & (letter != 'T'): 
            window = ''
        else:
            window = window+ letter
            counter = counter + 1
            if counter == klen:
                returnable +=  [ord_kmer(window)] #
                window = window[1:]
                counter = counter - 1
    return returnable


cdef ord_kmer_v3(str kmer, int val):
    cdef str to_ret = ''
    cdef int i = 0
    while i < val:
        inv = char_transform(kmer[val - i-1]) 
        if kmer[i]< inv:
            return kmer
        elif kmer[i] > inv:
            while i < val:
                to_ret = str(char_transform(kmer[val - i-1])) + to_ret
                i += 1
            return to_ret
        to_ret = str(i) + to_ret
        i +=1
    return kmer

def sklearn_ngram(text, n):
    cdef int text_len = len(text)
    cdef list ngrams = []
    ngrams_append = ngrams.append
    for i in range(text_len - n + 1):
        ngrams_append(ord_kmer_v2(text[i: i + n], n))
    return ngrams

def sklearn_ngram_yield(text, n):
    cdef int text_len = len(text)
    #cdef list ngrams = []
    for i in range(text_len - n + 1):
        yield ord_kmer_v2(text[i: i + n], n)
        
def sklearn_ngram_yield_native(text, n):
    cdef int text_len = len(text)
    #cdef list ngrams = []
    for i in range(text_len - n):
        yield ord_kmer_v2(text[i: i + n], n)

def ngram_native( filename,  klen = 20):
    forward  =  ''
    counter = 0
    sort_int = 0
    with open(filename, 'r') as f:
        while 1:
            if not f.readline(): break
            while 1:
                c = f.read(1) 
                if c == '\n' : break
                forward += c
                counter += 1
                if counter == klen:
                    yield ord_kmer_v2(forward,klen)
                    forward = forward[1:]
                    counter -= 1
            f.readline()
            f.readline()



def ngram_native_skl( filename, klen = Kmer_length_func):
    with open(filename, 'r') as f:
        while 1:
            if not f.readline(): break
            for i in sklearn_ngram_yield_native(f.readline(), klen): yield i
            f.readline()
            f.readline()   
            
def ngram_nltk( filename, klen = Kmer_length_func):
    with open(filename, 'r') as f:
        while 1:
            if not f.readline(): break
            for i in ngrams(f.readline(), klen): 
                yield ord_kmer_tup_py(i,klen)
            f.readline()
            f.readline()       

def Daskfilename2kmers(filename, blocksize= '1MiB', reversible = False):
    if not reversible: return db.read_text(filename, blocksize= blocksize).map(tokenize_read_NOORD).flatten().compute()
    return db.read_text(filename, blocksize= blocksize).map(tokenize_read).flatten().compute() #, blocksize='.5MiB'

def filename2kmers_pandas(filename, Klen = Kmer_length_func):
    return chain.from_iterable(pd.read_csv(filename, header = None, engine = 'c', dtype = str, names = ['read'])[1::4].read.apply(sklearn_ngram, args=(Kmer_length_func,)))

def filename2kmers_dd(filename, blocksize= '1MiB', Klen = Kmer_length_func):
    return chain.from_iterable(dd.read_csv(filename, header = None, engine = 'c', dtype = str, names = ['read'], blocksize= blocksize).read.apply(tokenize_read, args=(Kmer_length_func,),meta=('read', 'object')).compute())


def filename2kmers_dd_skl(filename):
    df = dd.read_csv(filename, header = None, engine = 'c', dtype = str, names = ['read'], blocksize= '1MiB').read.apply(sklearn_ngram, args=(Kmer_length_func,),meta=('read', 'object')).compute()[1::4]
    ret = list(chain.from_iterable(df))
    client.cancel(df)
    return ret

def filename2kmers_dd_skl_yield_yield(filename):
    return chain.from_iterable(dd.read_csv(filename, header = None, engine = 'c', dtype = str, names = ['read'], blocksize= '1MiB').read.apply(sklearn_ngram_yield, args=(Kmer_length_func,),meta=('read', 'object')).compute()[1::4])

def filename2kmers_explode(filename, blocksize= '1MiB'):
    return list(dd.read_csv(filename, blocksize= blocksize, header = None, engine = 'c', dtype = str, names = ['read']).read.apply(tokenize_read,meta=('read', 'object')).explode().compute())

def Daskfilenames_to_kmers_sparse(files = [], n_features=int(1e6), tfidf_bool = True, hasher = True, minmaxfreqkmer = (0.08,0.92), feature_names = False):
    if hasher:
        ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = Daskfilename2kmers).fit_transform(files)
        if tfidf_bool: return TfidfTransformer().fit_transform(ret)
        return ret
    if not hasher:
        vec = TfidfVectorizer(analyzer=Daskfilename2kmers, max_df = minmaxfreqkmer[1], min_df = minmaxfreqkmer[0], max_features = n_features )
        X = vec.fit_transform(files)
        if feature_names: return (X, np.array(vec.get_feature_names()))
        return X

    
def iterator_tokenize_read(reads):
    output = []
    for read in reads: output.append(tokenize_read(read))
    return output   

    
    
def Genomes2Sparse(files = [], n_features=int(1e6), tfidf = True, hasher = True, MAF = 0, feature_names = False, method = 'dask', wait = False, kmer_length = 20):
    global Kmer_length_func
    Kmer_length_func = kmer_length
    if hasher:
        if method == 'dask':
            ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = filename2kmers_dd_skl_yield_yield).fit_transform(db.from_sequence(files)).compute()
        if method == 'pandas':
            ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = filename2kmers_pandas).fit_transform(db.from_sequence(files)).compute()
        if method == 'native':
            if wait: ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = ngram_native_skl).fit_transform(db.from_sequence(files)).compute()
            else:
                tokens = client.map(ngram_native_skl, files)
                ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = yielder_).fit_transform(client.gather(tokens))
            
        if MAF > 0:
            minimum,maximum = sorted([MAF, 1-MAF])
            n_doc = len(files)
            ret = ret[:, (minimum*n_doc <ret.getnnz(0)) & (ret.getnnz(0)< (maximum)*n_doc )]
            if ret.shape[1] == 0: return 'no Kmer with MAF of interest'
        if tfidf: ret =  TfidfTransformer().fit_transform(ret)
        return ret
    
    if not hasher:
        minimum,maximum = sorted([MAF, 1-MAF])
        if dask: vec = TfidfVectorizer(analyzer=filename2kmers_dd_skl_yield_yield, max_df = maximum, min_df = minimum, max_features = n_features )
        if not dask:  vec = TfidfVectorizer(analyzer=filename2kmers_pandas, max_df = maximum, min_df = minimum, max_features = n_features )
        X = vec.fit_transform(files)
        if feature_names: return (X, np.array(vec.get_feature_names()))
        return X
    
def Genomes2Sparse_delayed(files = [], n_features=int(1e6), tfidf = True, hasher = True, MAF = 0, feature_names = False, kmer_length = 20):
    global Kmer_length_func
    Kmer_length_func = kmer_length
    tokens = client.map(filename2kmers_dd, files)
    if hasher:
        ret = dask_text.HashingVectorizer(n_features=n_features, analyzer = lambda x: x).fit_transform(client.gather(tokens))
        if MAF > 0:
            minimum,maximum = sorted([MAF, 1-MAF])
            n_doc = len(files)
            ret = ret[:, (minimum*n_doc <ret.getnnz(0)) & (ret.getnnz(0)< (maximum)*n_doc )]
            if ret.shape[1] == 0: return 'no Kmer with MAF of interest'
        if tfidf: ret =  TfidfTransformer().fit_transform(ret)
        return ret
    if not hasher:
        minimum,maximum = sorted([MAF, 1-MAF])
        vec = TfidfVectorizer(analyzer=lambda x:x, max_df = maximum, min_df = minimum, max_features = n_features )
        X = vec.fit_transform(client.gather(tokens))
        if feature_names: return (X, np.array(vec.get_feature_names()))
        return X
    


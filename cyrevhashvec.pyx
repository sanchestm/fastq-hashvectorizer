from libc.stdlib cimport malloc, free
from libc.stdio cimport fopen, fclose, FILE, EOF, fgetc, feof
from scipy.sparse import coo_matrix
import numpy as np

cdef unsigned char char_rev(unsigned char chars):
    if chars == 65: return 84
    elif chars == 67: return 71
    elif chars == 71: return 67
    elif chars == 84: return 65
    return 78

cdef unsigned char char2int(char chars):
    if chars == 65: return 0
    elif chars == 67: return 1
    elif chars == 71: return 2
    elif chars == 84: return 3
    return -1

cdef char check_sort2(unsigned char * fw, unsigned char * rv, int klen):
    cdef int i
    for i in range(klen):
        if fw[i] < rv[i]: return 1
        elif fw[i] > rv[ i]: return 0
    return 1

cdef int RevHash(unsigned char *kmer):
    cdef unsigned char c
    cdef int rev = 0
    cdef int i
    for i in range(len(kmer)):
        rev = rev*4 + char2int(kmer[i])
    return rev

cdef long longRevHash(unsigned char *kmer):
    cdef unsigned char c
    cdef long rev = 0
    cdef int i
    for i in range(len(kmer)):
        rev = rev*4 + char2int(kmer[i])
    return rev

def longRevHashPY(kmer):
    return longRevHash(kmer)

cdef list cycols(str filename, int klen):
    """Efficiently read in a file"""
    cdef FILE *fp = NULL # create a file pointer
    fp = fopen(filename.encode(encoding='utf-8'), "rb")
    if fp == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)
    # file parsing variables
    cdef char c 
    cdef int i
    cdef int is_line = 1 
    cdef int window_pos = 0
    cdef bytes pywindow 
    cdef unsigned char *cwindow = <unsigned char *> malloc((klen+1)*sizeof(char))
    cdef unsigned char *invcwindow = <unsigned char *> malloc((klen+1)*sizeof(char))
    cdef list col = []

    #cdef int stopearly = 0
    # bypass the gil
    #with nogil:
    while 1 :
        c = fgetc(fp)
        if feof(fp): break
        if c == 78: window_pos = 0
        elif c == 10: 
            is_line += 1
            if is_line > 4: is_line= 1
            window_pos = 0
        elif is_line == 2: 
            cwindow[window_pos] = c
            invcwindow[klen - window_pos-1] = char_rev(c)
            window_pos += 1
            if window_pos == klen: 
                if check_sort2(cwindow, invcwindow, klen ) == 1: 
                    col.append(longRevHash(cwindow))
                else: col.append(longRevHash(invcwindow))
                for i in range(klen-1): cwindow[i] = cwindow[i+1]
                for i in range(klen-2,-1,-1): invcwindow[i+1] = invcwindow[i]
                window_pos -= 1
    # close the file
    free(cwindow)
    free(invcwindow)
    fclose(fp)
    return col

def cyhashvectorizer(list filenames = [], int klen = 20, memorysafe = True):
    cdef long size = 4**klen - 4**(klen//2)  #cdef long long 
    cdef int corpus_size = len(filenames)
    ret = coo_matrix(([], ([],[])),shape=(corpus_size, size), dtype = 'int')
    for  num,name  in enumerate(filenames):
        col = cycols(name, klen)
        length = len(col)
        ret += coo_matrix((np.ones(length, dtype = int), (np.full(length, num, dtype=int), col)), shape=(corpus_size, size), dtype = 'int')
        if memorysafe: ret.sum_duplicates()
    return ret.tocsr()



import numpy as np   
import scipy.sparse as sp

from libc.stdlib cimport malloc, free
from libc.stdio cimport fopen, fclose, FILE, EOF, fgetc, feof
from cython.parallel import prange
cimport numpy as cnp
cimport cython

cdef int char2int(char chars) nogil:
    if chars == 65: return 0
    elif chars == 67: return 1
    elif chars == 71: return 2
    elif chars == 84: return 3
    return -1

cdef int invchar2int(char chars) nogil:
    if chars == 65: return 3
    elif chars == 67: return 2
    elif chars == 71: return 1
    elif chars == 84: return 0
    return -1


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)
cdef cykmerarray(str filename, int klen, unsigned long max_size = -1): ##unsigned long is up to 32bp cnp.uint16_t[:]
    """Efficiently read in a file"""
    cdef FILE *fp = NULL # create a file pointer
    fp = fopen(filename.encode(encoding='utf-8'), "rb")
    if fp == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)
    # file parsing variables
    cdef unsigned long array_size = 4**klen - 4**(klen//2)
    cdef unsigned long fw_div = 4**(klen-1)
    
    if max_size > 0 and max_size < array_size:
        array_size = max_size
    
    cdef char c 
    cdef int i
    cdef int is_line = 1 
    cdef int window_pos = 0
    cdef bytes pywindow 
    cdef unsigned long fw_num = 0
    cdef unsigned long rv_num = 0
    #cdef list col = []
    cdef cnp.ndarray[cnp.int16_t, ndim=1] result = np.zeros((array_size), dtype = np.int16)
    
    #cdef int stopearly = 0
    # bypass the gil
    with nogil:
        while 1 :
            c = fgetc(fp)
            if feof(fp): break
            if c == 78: window_pos = 0
            elif c == 10: 
                is_line += 1
                if is_line > 4: is_line= 1
                window_pos = 0
            elif is_line == 2: 
                fw_num = fw_num*4 + char2int(c)
                rv_num = invchar2int(c)*(4**window_pos) + rv_num
                window_pos += 1
                if window_pos == klen: 
                    if fw_num < rv_num: result[fw_num%max_size] += 1
                    else: result[rv_num%max_size] += 1
                    fw_num %= fw_div
                    rv_num //= 4
                    window_pos -= 1
    # close the file
    fclose(fp)
    return result

cdef list getK(char *string, int klen):
    cdef unsigned long fw_num = 0
    cdef unsigned long rv_num = 0
    cdef unsigned long fw_div = 4**(klen-1)
    cdef window_pos = 0
    cdef char c
    ret = []
    for c in string:
        if char2int(c) == -1: window_pos = 0
        else:
            fw_num = fw_num*4 + char2int(c)
            rv_num = invchar2int(c)*(4**window_pos) + rv_num
            window_pos += 1
            if window_pos == klen: 
                    if fw_num < rv_num: ret.append(fw_num)
                    else: ret.append(rv_num)
                    fw_num %= fw_div
                    rv_num //= 4
                    window_pos -= 1
    return ret

def pygetK(string, klen):
    return getK(string.encode(encoding='ascii'), klen)
    

def cyhashSMARTVEC(filenames = [], int klen = 13, max_size = -1):
    if max_size > 1e9: 
        print('this will use more than 1 gb per file')
        return -1
    cdef long size = 4**klen - 4**(klen//2)  #cdef long long 
    if max_size > 0 and max_size < size:
        size = max_size
    cdef int corpus_size = len(filenames)
    for  num,name  in enumerate(filenames):
        if num == 0: ret = sp.csr_matrix(cykmerarray(name, klen, size))
        else: ret = sp.vstack((ret, cykmerarray(name, klen, size)))
    ret._meta = sp.eye(0, dtype = np.int16 ,format="csr")
    return ret
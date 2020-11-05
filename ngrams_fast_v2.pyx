from libc.stdlib cimport malloc, free
from libc.stdio cimport fopen, fclose, FILE, EOF, fgetc, feof

cdef char char_trans(char chars):
    if chars == 65: return 84
    elif chars == 67: return 71
    elif chars == 71: return 67
    elif chars == 84: return 65
    return 78

cdef char check_sort(char * fw, char * rv, int klen):
    cdef int i
    for i in range(klen):
        if fw[i] < rv[i]: return 1
        elif fw[i] > rv[i]: return 0
    return 1

def cyngrams(str filename, int klen):
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
    cdef  char *cwindow = < char *> malloc((klen+1)*sizeof(char))
    cdef  char *invcwindow = < char *> malloc((klen+1)*sizeof(char))
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
            invcwindow[klen - window_pos-1] = char_trans(c)
            window_pos += 1
            if window_pos == klen: 
                if check_sort(cwindow, invcwindow, klen ) == 1: yield cwindow#.decode(encoding="ascii",errors="ignore")
                else: yield invcwindow#.decode(encoding="ascii",errors="ignore")
                #yield cwindow
                #yield invcwindow
                for i in range(klen-1): cwindow[i] = cwindow[i+1]
                for i in range(klen-2,-1,-1): invcwindow[i+1] = invcwindow[i]
                window_pos -= 1
    # close the file
    free(cwindow)
    free(invcwindow)
    fclose(fp)
    return

#def Cyngrams(filename, klen = 20):
#    for i in cyngrams(filename, klen): yield (i)
        
def listcyngrams(filename, klen = 20):
    return list(cyngrams(filename, klen))

cdef list listcyngrams_fast(str filename, int klen):
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
    cdef char *cwindow = < char *> malloc((klen+1)*sizeof(char))
    cdef char *invcwindow = < char *> malloc((klen+1)*sizeof(char))
    cdef list ret = []
    
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
            invcwindow[klen - window_pos-1] = char_trans(c)
            window_pos += 1
            if window_pos == klen: 
                if check_sort(cwindow, invcwindow, klen ) == 1: ret.append(cwindow)#.decode(encoding="ascii",errors="ignore")
                else: ret.append(invcwindow)#.decode(encoding="ascii",errors="ignore")
                #yield cwindow
                #yield invcwindow
                for i in range(klen-1): cwindow[i] = cwindow[i+1]
                for i in range(klen-2,-1,-1): invcwindow[i+1] = invcwindow[i]
                window_pos -= 1
    # close the file
    free(cwindow)
    free(invcwindow)
    fclose(fp)
    return ret


def listcyngrams_fastpy(filename, klen = 20):
    return listcyngrams_fast(filename, klen)
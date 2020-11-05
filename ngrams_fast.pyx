tab_b = bytes.maketrans(b"ACGT", b"TGCA")

def bRevComp(seq):
    return seq.translate(tab_b)[::-1]

cdef int char_transform(int char):
    if char == 65: return 84
    elif char == 67: return 71
    elif char == 71: return 67
    elif char == 84: return 65
    return 78

cdef OK(unsigned char *kmer, int val):
    cdef int i = 0
    while i < val:
        if kmer[i] <char_transform(kmer[val - i-1]):
            return kmer
        elif kmer[i] > char_transform(kmer[val - i-1]):
            return bRevComp(kmer)
        i +=1
    return kmer


def sklearn_ngram_yield_fast(unsigned char *text, int n):
    cdef int text_len = len(text) 
    cdef int i = 0 
    for i in range(text_len - n):
        yield OK(text[i: i + n], n)
        
def sklearn_ngram_yield_fast_old(text, int n):
    cdef int text_len = len(text)
    for i in range(text_len - n):
        yield OK(text[i: i + n], n)

def yieldkmers( filename, klen = 20):
    with open(filename, 'rb') as f:
        while 1:
            if not f.readline(): break
            for i in sklearn_ngram_yield_fast_old(f.readline(), klen): yield i
            f.readline()
            f.readline()  
            
def listkmers(filename, klen = 20):
    return list(yieldkmers(filename, klen))
            

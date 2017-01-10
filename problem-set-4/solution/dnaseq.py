#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        self.dict = {}
        for pair in pairs:
            self.put(pair[0], pair[1])

    # Associates the value v with the key k.
    def put(self, k, v):
        if self.dict.has_key(k):
            self.dict[k].append(v)
        else:
            self.dict[k] = [v]

    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        return self.dict[k] if self.dict.has_key(k) else []

# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    subseq = ''
    i = 0
    rh = None
    hash_init = False
    previtm = ''
    for c in seq:
        subseq = ''.join([subseq, c])
        if len(subseq) % k == 0:
            if not hash_init:
                rh = RollingHash(subseq)
                hash_init = True
            else:
                rh.slide(previtm, c)   
            yield (rh.current_hash(), (i-k+1, subseq))
            previtm = subseq[0]
            subseq = subseq[1:] #Get rid of the first char
        i += 1

# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    subseq = ''
    i = 0
    rh = None
    hash_init = False
    previtm = ''
    for c in seq:
        subseq = ''.join([subseq, c])
        if len(subseq) % k == 0:
            if not hash_init:
                rh = RollingHash(subseq)
                hash_init = True
            else:
                rh.slide(previtm, c)   
            yield (rh.current_hash(), (i-k+1, subseq))
            previtm = subseq[0]
            subseq = subseq[1:] #Get rid of the first char
            #skip the next m-k elements
            for j in range(0, m-k):
                seq.next()
            i += m-k
        i += 1

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    multidict_a = Multidict(intervalSubsequenceHashes(a, k, m))

    for hash_b, (b_offset, b_subseq) in subsequenceHashes(b, k):
        for (a_offset, a_subseq) in multidict_a.get(hash_b):
            if a_subseq == b_subseq:
                yield (a_offset, b_offset)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0])
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)

import random

# Hash function that uses xor properties (from https://github.com/yoann-dufresne/bloomtest/blob/master/hash.cpp)
# input: x -> value to hash
# output: hashed value
def xorshift64(x):
    x ^= x << 13
    x ^= x >> 7
    x ^= x << 17
    return x

# Generate multiple hashes on a value in using the xorshift function
# input: bf -> object of BloomFilter class
#        x  -> value to hash
#        k  -> size of kmers
# output: index generator indicating a cell of the Bloom filter and the corresponding bit
def multihash(bf, x, k):
    if x==0:
        x=sum([2**i for i in range(1, 2*k, 2)])  # replace x by the code of the reverse sequence if x is 0
    for i in range(bf.nf):
        x = xorshift64(x)
        yield x%bf.size, x%8

# Class of Bloom Filter
# attributes: n    ->  size of the filter
#             nf   ->  number of hashing function to apply
#             size ->  size of the array
#             array  ->  array containing the filter (a binary value in each cell)
class BloomFiltrer:
    
    # Constructor
    # input: n  -> size of the filter
    #        nf -> number of hashing function to apply
    def __init__(self, n, nf):
        self.n = n
        self.nf = nf
        self.size = int(n/8)+1
        self.array = bytearray(self.size)
    
    # Add a value to the filter
    # input: x -> value to add
    #        k -> size of kmers
    def add_value(self, x, k):
        for ind,bit in multihash(self, x, k):
            self.array[ind]|=1<<bit
    
    # Verify if a value is in the filter
    # input: x -> value to add
    #        k -> size of kmers
    # output: True if x is in the filter, False otherwise
    def is_present(self, x, k):
        for ind,bit in multihash(self, x, k):
            if self.array[ind]&(1<<bit)==0: return False
        return True
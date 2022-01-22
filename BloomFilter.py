from src.Bloom_Filter_class import *
from src.Hashing_sequence import *
import sys, os


if len(sys.argv)!=6:
    raise ValueError('Not enough parameters')
fasta = sys.argv[1]     # path to the fasta file
k = int(sys.argv[2])    # size of kmers
if k<1 or k>31:
    raise ValueError('k not between 1 and 31')
n = int(sys.argv[3])    # size of the filter
nf = int(sys.argv[4])   # number of hash function to apply
if nf<1:
    raise ValueError('nf not strictly positive')
r = int(sys.argv[5])    # number of requests to test on the filter
if r<0:
    raise ValueError('r not positive')


# Initialization
bf = BloomFiltrer(n, nf)
seq = iter(read_fasta(fasta))
code, reverse_code = hashing(''.join([next(seq) for i in range(k)]))
bf.add_value(min(code, reverse_code), k)

# Iteration of Bloom Filter construction
for nt in seq:
    code, reverse_code = next_kmer(code, reverse_code, nt, k)
    bf.add_value(min(code, reverse_code), k)

# Test the random requests
cpt = 0
for i in range(r):
    if bf.is_present(random.randrange(2**(2*k-1)), k):
        cpt+=1
print(cpt)
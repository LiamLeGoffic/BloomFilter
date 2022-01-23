from src.Bloom_Filter_class import *
from src.Hashing_sequence import *
import sys


if len(sys.argv)!=6:
    raise ValueError('6 parameters required')
fasta = sys.argv[1]     # path to the fasta file
if fasta[-6:]!='.fasta':
    raise ValueError('The File does not have the fasta format')
k = int(sys.argv[2])    # size of kmers
if k<1 or k>31:
    raise ValueError('k not between 1 and 31')
n = int(sys.argv[3])    # size of the filter
if n > 17179869184 or n<1:
    raise ValueError('n not between 1 and 17179869184')
nf = int(sys.argv[4])   # number of hash function to apply
if nf<1 or nf>64:
    raise ValueError('nf between 1 and 64')
r = int(sys.argv[5])    # number of requests to test on the filter
if r<0:
    raise ValueError('r not positive')


# Initialization
bf = BloomFiltrer(n, nf)
seq = iter(read_fasta(fasta))
code, reverse_code = hashing(''.join([next(seq) for i in range(k)]))
bf.add_value(min(code, reverse_code), k)

# Iterations of the Bloom Filter construction
for nt in seq:
    code, reverse_code = next_kmer(code, reverse_code, nt, k)
    bf.add_value(min(code, reverse_code), k)

# Test the random requests
cpt = 0
for i in range(r):
    if bf.is_present(random.randrange(4**k), k):
        cpt+=1
print(cpt)

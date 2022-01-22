# Read a fasta file nucleotide by nucleotide
# input:  file ->  path to fasta file
# output: generator of nucleotide in string format
def read_fasta(file):
    fasta = open(file, 'r')
    for line in fasta:
        if line[0]=='>':
            continue
        for nt in line:
            if nt in ['N', '\n']:
                continue
            yield nt

# Hash a sequence of nucleotides and its reverse sequence by replacing A, C, T and G with 0, 1, 2 and 3 respectively
# input:  seq  ->  nucleotide sequence in string format
# output: hashed sequence and hashed reverse sequence
def hashing(seq):
    code, r_code = 0, 0
    n = len(seq)
    for i in range(n):
        nt = (ord(seq[i])>>1)&0b11
        comp = ((~nt)&0b10)+(nt&0b01)
        code = (code<<2)+nt
        r_code = ((r_code>>2)+(comp<<((n-1)*2)))
    return code, r_code

# Get the code of the next kmer thanks to the last one and the next nucleotide
# input: kmer      ->  last hashed kmer
#        r_kmer    ->  last hashed reverse kmer
#        next_kmer ->  next nucleotide character
#        k         ->  size of kmers
# output:  next hashed kmer code and next hashed reverse kmer
def next_kmer(kmer, r_kmer, next_nt, k):
    next_nt = (ord(next_nt)>>1)&0b11
    comp = ((~next_nt)&0b10)+(next_nt&0b01)
    mask = (1<<(2*(k-1)))-1
    return ((kmer&mask)<<2)+next_nt, (r_kmer>>2)+(comp<<((k-1)*2))
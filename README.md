# Bloom Filter
## What is it ?
Wikipedia page : https://fr.wikipedia.org/wiki/Filtre_de_Bloom

## What does the algorithm do ?
The algorithm builds a Bloom filter from the kmers of a nucleotide sequence (from a FASTA file) and test a defined number of requests (corresponding to random kmers) on the filter to show if they belong to the sequence. The number of requests that pass the filter is returned.

For more informations about the project instructions : https://github.com/yoann-dufresne/bloomtest

## Usage
Five parameters are required :
- path : path that leads to the FASTA file
- k : size of kmers (max 31)
- n : size of the filter (max 2^34 (16 Go))
- nf : number of hashing function (max 64)
- r : number of requests

To run, use this command :
```
python BloomFilter.py -path -k -n -nf -r
```

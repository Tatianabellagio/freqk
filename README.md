# freqk

Estimate frequencies of known variants in pool-seq data from k-mer counts

## Inputs

1. Reference path in fasta format

2. Index of fasta file `samtools faidx <ref fasta>`

3. bgzipped VCF file of variants called against reference, variants should be phased and non-overlapping

4. Index of vcf file

## Usage

1. Index the panel of reference variants

`freqk index --vcf tests/1.vcf.gz --fasta tests/1.fasta --k 31 --output index.txt`

2. Count the k-mers specified in the index

3. Convert k-mer counts to allele frequencies

## To-do

- [ ] for count subcommand: load in k-mers from index as hash set, loop over input reads with a fastx iterator, slide window to get k-mers, only count k-mers if they're in the hash set

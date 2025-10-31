# freqk

Estimate frequencies of known variants in pool-seq data from k-mer counts

## Inputs

1. Reference path in fasta format

2. Index of fasta file `samtools faidx <ref fasta>`

3. bgzipped VCF file of variants called against reference, variants should be phased and non-overlapping

4. Index of vcf file

## Usage

1. Index the panel of reference variants

`freqk index --vcf tests/1.vcf.gz --fasta tests/1.fasta -k 31 --output index.txt`

2. Count the k-mers specified in the index

`freqk count --index index.txt --reads tests/1_R1.fastq tests/1_R2.fastq --output counts.txt`

3. Call allele frequences based on k-mer counts

`freqk call --index index.txt --counts counts.txt --output calls.txt`

## To-do

- [x] parallelize for loop over reads

- [x] remove leftover code

- [x] organize code into modules

- [x] use paraseq instead of fastqrs - paraseq has a conflict that prevents it from working with rusthtslib, also seems to do similarly to fastqrs in benchmarks

- [ ] convert counts into allele frequencies

- [ ] deduplicate: remove any allele-specific k-mers found elsewhere in the reference genome

- [x] dedeuplicate index: remove any k-mers that are shared across variants

- [ ] add filter or warning for when 0 or only a few k-mers tag a variant (if not many k-mers tag a variant and the variant is rare in the pool, then coverage needs to be absurdly high in order to accurately estimate allele frequency)

- [ ] add unit tests

- [x] determine k-mer length from index

- [x] when indexing, count number of k-mers for each allele

- [x] for count subcommand: load in k-mers from index as hash set, loop over input reads with a fastx iterator, slide window to get k-mers, only count k-mers if they're in the hash set

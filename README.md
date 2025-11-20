# freqk

Estimate frequencies of known variants in pool-seq data from k-mer counts

## Installation

Grab rust binary from release page

## Inputs

1. Reference sequence in fasta format

2. Index of fasta file `samtools faidx <ref fasta>`

3. VCF file of variants called against reference, should be 

* sorted

* bgzipped

* only non-overlapping

* genome-wide

4. Index of vcf file

5. Pooled DNA sequencing reads

## Usage

1. Index the panel of reference variants

`freqk index -f tests/1.fasta -v tests/1.vcf.gz -o index.txt -k 31`

2. (Optional but recommended) Deduplicate index

This step can be done in 2 stages. You can do one, both, or neither stages, but generally doing both before proceeding to the counting step (3) is strongly recommended for best results.

`freqk var-dedup --index index.txt --output dedup.txt`

Another round of deduplication will remove any putatively allele-specific k-mers that are found elsewhere in the reference sequence. 

`freqk ref-dedup -i var_index.txt -o ref_index.txt -f tests/1.fasta -v tests/1.vcf.gz`

3. Count the indexed k-mers in the pool-seq reads using two threads (`-n 2`). 

`freqk count -i ref_index.txt -r tests/all.fastq.gz -n 2 -f count_by_allele.txt -c count_by_kmer.txt`

4. Normalize allele-specific k-mer counts into allele frequencies

This step just divides the counts of allele-specific k-mers in the reads by the number of allele-specific k-mers in the index.

`freqk call -i ref_index.txt -c counts_by_allele.txt -o calls.txt`

## To-do

- [x] get k-mer length from index in ref-dedup

- [x] check that REF allele matches fasta file

- [ ] add typical counting speed (number of 150 bp reads per second per thread) to README

- [ ] add input checkers, check that vcf is sorted for freqk index

- [ ] add unit tests

- [ ] add more methods to structs

- [ ] add filter or warning for when 0 or only a few k-mers tag a variant (if not many k-mers tag a variant and the variant is rare in the pool, then coverage needs to be absurdly high in order to accurately estimate allele frequency)

- [x] parallelize for loop over reads

- [x] remove leftover code

- [x] organize code into modules

- [x] use paraseq instead of fastqrs - paraseq has a conflict that prevents it from working with rusthtslib, also seems to do similarly to fastqrs in benchmarks

- [x] convert counts into allele frequencies

- [x] skip over variants at chromosome ends and overlapping variants 

- [x] check that reference sequence in vcf matches fasta

- [x] var-dedup: add more print messages

- [x] bug: var-dedup, variant at end of one chromosome and start of another chromosome are labeled as overlapping and skipped

- [x] bug: ref-dedup, lots of variants are skipped

- [x] deduplicate: remove any allele-specific k-mers found elsewhere in the reference genome

- [x] dedeuplicate index: remove any k-mers that are shared across variants

- [x] determine k-mer length from index

- [x] when indexing, count number of k-mers for each allele

- [x] for count subcommand: load in k-mers from index as hash set, loop over input reads with a fastx iterator, slide window to get k-mers, only count k-mers if they're in the hash set

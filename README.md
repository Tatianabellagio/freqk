# freqk

Estimate frequencies of known variants in pool-seq data from k-mer counts

Author: Miles Roberts

Contact: leave a github issue

**Note:** Before leaving an issue, re-run the problematic command(s) with `-vvv` to get full debug information for that command and post the debug information with your issue.

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

## Outputs

Estimates of allele frequencies for variants in the vcf file. These are output in a text file with pipe separators

```
0.01|0.99
1|0
0.75|0.25
0.3|0.5|0.2
```

The number of alleles per site varies, so the number of entries per line varies. If you want to load this file in R. You can do something like the following

## Quick start

1. Index the panel of reference variants

`freqk index -f tests/1.fasta --vcf tests/1.vcf.gz -o index.txt -k 31`

2. (Optional but recommended) Deduplicate index

This step can be done in 2 stages. You can do one, both, or neither stages in any order, but generally doing both before proceeding to the counting step (3) is strongly recommended for the most rigorous results.

The faster of the two deduplication steps is `var-dedup`. This simply scans through the index and removes any allele-specific k-mers that are found at other variants in the index. For example:

`freqk var-dedup --index index.txt --output dedup.txt`

The slower deduplication step is `ref-dedup`. This step removes any allele-specific k-mers that are found elsewhere in the reference sequence. 

`freqk ref-dedup -i var_index.txt -o ref_index.txt -f tests/1.fasta --vcf tests/1.vcf.gz`

3. Count the indexed k-mers in the pool-seq reads (multithreaded). 

For example, counting indexed k-mers with four threads (`-n 4`) looks like this:

`freqk count -i ref_index.txt -r tests/all.fastq.gz -n 4 -f count_by_allele.txt -c count_by_kmer.txt`

4. Normalize allele-specific k-mer counts into allele frequencies

This step just divides the counts of allele-specific k-mers in the reads by the number of allele-specific k-mers in the index.

`freqk call -i ref_index.txt -c counts_by_allele.txt -o calls.txt`

## Command line interfaces

All commands have a verbosity flag. Only errors are output by default, but adding `-v` will make warnings print, `-vv` means info will also print, and `-vvv` means debug data will print.

## To-do

- [ ] hash info for variants so that you skip duplicates?

- [ ] for ref-dedup, what to do if there are chromosomes in fasta that have no variants in vcf? 

- [ ] for ref-dedup, put k-mers into hashset at end of each loop?

- [x] add verbose flag with clap, improve logging: https://rust-cli.github.io/book/tutorial/output.html

- [ ] index step, for variants within k bp, only get k-mers that overlap one variant

- [x] skip invariant sites, sites where REF and ALT are the same or there is no ALT information

- [x] get k-mer length from index in ref-dedup

- [x] check that REF allele matches fasta file

- [ ] add typical counting speed (number of 150 bp reads per second per thread) to README

- [ ] add input checkers, check that vcf is sorted for freqk index

- [ ] add unit tests

- [ ] add more methods to structs

- [x] add filter or warning for when 0 or only a few k-mers tag a variant (if not many k-mers tag a variant and the variant is rare in the pool, then coverage needs to be absurdly high in order to accurately estimate allele frequency)

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

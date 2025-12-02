use rust_htslib::bcf::Reader;
use rust_htslib::bcf::Read;
use bio::io::fasta::IndexedReader;
use std::fs::File;
use std::io::{Write, BufWriter};
use bio::bio_types::genome::AbstractLocus;
use std::collections::HashSet;
use std::collections::HashMap;
use crate::common;

// given a list of k-mers by allele
// find all k-mers shared across alleles
fn find_dup_kmers(mut data: Vec<Vec<String>>) -> Vec<Vec<String>> {
    // first pass: count how many times k-mers are found across alleles
    let mut counts: HashMap<String, usize> = HashMap::new();

    for inner_vec in &data {
        for s in inner_vec {
            *counts.entry(s.to_string()).or_insert(0) += 1;
        }
    }

    // identify k-mers found more than once
    let dup_kmers: Vec<String> = counts
        .into_iter() 
        .filter(|(_key, value)| *value > 1) // Filter pairs where k-mer is non-unique
        .map(|(key, _value)| key) // get only the keys
        .collect(); // Collect the keys into a Vec

    let dup_kmers_hashset: HashSet<String> = dup_kmers.into_iter().collect();

    // second pass: remove any k-mers that were duplicated
    for inner_vec in &mut data {
        inner_vec.retain(|s| !dup_kmers_hashset.contains(s));
    }

    data
}


// insert variant into reference, get k-mers for each variant
pub fn index_workflow(vcf_path: &String, fasta_path: &String, output_path: &String, k: &i64) -> Option<Vec<i32>> {
    // rust-htslib provides VCF I/O.
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");

    // read indexed fasta file
    let mut faidx = IndexedReader::from_file(fasta_path).unwrap();

    // read fasta index to get chrom_lengths
    let chrom_lengths = common::read_fai(fasta_path);
    log::info!("Chromosome lengths:");
    log::info!("{:?}", chrom_lengths);

    // output file
    let mut buffered_file = BufWriter::new(File::create(output_path).ok()?);

    // iterate through each row of the vcf body.i
    let mut i = 0;
    let mut vcf_iterator = vcf_reader.records().peekable();
    while let Some(record_result) = vcf_iterator.next() {
        i += 1;
        let record = record_result.expect("Failed to read record!");
        let pos = record.pos() - 1;
        let chrom = record.contig();
        let mut ku = *k as usize;
        log::debug!("Processing record CHROM: {} POS: {}", chrom, pos);
        // get putative region, but check edge cases
        let mut region_start = pos - k + 1;
        let region_end = pos + k;
        // check if variant near chromosome ends or not found in reference
        if let Some(num_ref) = chrom_lengths.as_ref().expect("Failed to read .fa.fai file!").get(chrom) {
            let end = *num_ref;
            if end < *k {
                log::warn!("Chromosome shorter than k, skipping CHROM: {} POS: {}", chrom, pos);
                continue
            }
            if pos >= (end - k){
                log::warn!("Variant near chromosome end detected, skipping CHROM: {} POS: {}", chrom, pos);
                continue
                //region_end = end;
                //ku = 
            }
        } else {
            log::warn!("Warning: {} not found in .fa.fai file! Thus, will skip CHROM: {} POS: {}", chrom, chrom, pos);
            continue
        }
        if pos <= *k {
            log::warn!("Variant near chromosome start detected, skipping CHROM: {} POS: {}", chrom, pos);
            //continue
            region_start = 1;
            ku = (pos + 4) as usize;
        }
        // until you reach the end of the next file, peek at the next vcf record and see if it
        // overlaps with the current record, if there is an overlap then skip both this and the
        // next record
        if let Some(next_ref) = vcf_iterator.peek() {
            let next_result = next_ref.as_ref().unwrap();
            let pos_next = next_result.pos() - 1;
            let chrom_next = next_result.contig();
            if ((pos_next - pos) <= *k)  && (chrom == chrom_next) {
                log::warn!("Overlapping pair of variants detected, skipping CHROM: {} POS: {} and CHROM: {} POS: {}", chrom, pos, chrom_next, pos_next);
                vcf_iterator.next();
                continue
            }
        }
        // construct sequences for alleles
        log::debug!("Extracting allele sequences from VCF...");
        let mut alleles = String::new();
        for allele in record.alleles() {
            for c in allele {
                alleles.push(char::from(*c))
             }
            alleles.push(' ')
        }
        // Split the string by whitespace and collect into a Vec<&str>
        let alleles_list: Vec<&str> = alleles.split_whitespace().collect();
        // check if site is not variable
        //if alleles_list[0] == alleles_list[1] {
        if (alleles_list[1..alleles_list.len()]).contains(&alleles_list[0]) {
            log::warn!("Invariant site detected, skipping CHROM: {} POS: {}", chrom, pos);
            continue
        }
        // move the pointer in the index to the desired sequence and interval
        faidx.fetch(chrom, region_start.try_into().unwrap(), region_end.try_into().unwrap() ).expect("Couldn't fetch interval");
        // read the subsequence defined by the interval into a vector
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Couldn't read the interval");
        // convert to string
        let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
        // Split the string by whitespace and collect into a Vec<&str>
        //let alleles_list: Vec<&str> = alleles.split_whitespace().collect();
        // insert alleles into reference sequence to get variable sequences 
        // replace reference allele with variant allele
        let mut var_seqs = Vec::new();
        //let ku = *k as usize;
        let ref_allele_len = alleles_list[0].len();
        for allele in &alleles_list {
            let mut var_seq = seq_string.clone();
            var_seq.replace_range(ku..ku+ref_allele_len, allele);
            var_seqs.push(var_seq);
        }
        // check if REF allele matches fasta file
        if var_seqs[0] != seq_string {
            log::error!("REF allele does not match FASTA at CHROM: {} POS: {} \n FASTA: {} \n VCF  : {} \n REF: {}. Was the FASTA used as the reference to make the VCF?", chrom, pos, &var_seqs[0], &seq_string, &alleles_list[0]);
        }
        // get k-mers
        log::debug!("Extracting canonical k-mers from alleles...");
        let mut kmers_by_allele = Vec::new();
        for var_seq in &var_seqs {
            let kmer_list: Vec<String> = common::get_canonical_kmers(var_seq, *k as usize);
            kmers_by_allele.push(kmer_list);
        }
        // remove kmers shared across alleles
        log::debug!("Dedupping k-mers shared across alleles of the same locus...");
        let kmers_by_allele_no_dups = find_dup_kmers(kmers_by_allele);
        let mut joined_kmers_list: Vec<String> = Vec::new();
        let mut num_kmers_per_allele: Vec<String> = Vec::new();
        for inner_vec in kmers_by_allele_no_dups {
            let inner_joined = inner_vec.join(";");
            joined_kmers_list.push(inner_joined);
            num_kmers_per_allele.push(inner_vec.len().to_string());
        }
        // write index to output file
        let parts = vec![i.to_string(), chrom.to_string(), pos.to_string(), seq_string.clone(), alleles_list.join("|"), var_seqs.join("|"), num_kmers_per_allele.join("|"), joined_kmers_list.join("|")];
        writeln!(buffered_file, "{}", parts.join(",")).ok()?;
    }
    None
}


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
pub fn insert_var(vcf_path: &String, fasta_path: &String, output_path: &String, k: &i64) -> Option<Vec<i32>> {
    // rust-htslib provides VCF I/O.
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");

    // read indexed fasta file
    let mut faidx = IndexedReader::from_file(fasta_path).unwrap();

    // output file
    let mut buffered_file = BufWriter::new(File::create(output_path).ok()?);

    // iterate through each row of the vcf body.
    for (i, record_result) in vcf_reader.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut alleles = String::new();
        for allele in record.alleles() {
            for c in allele {
                alleles.push(char::from(*c))
             }
            alleles.push(' ')
        }

        // move the pointer in the index to the desired sequence and interval
        let pos = record.pos() - 1;
        let chrom = record.contig();
        faidx.fetch(chrom, (pos - k + 1).try_into().unwrap(), (pos + k).try_into().unwrap() ).expect("Couldn't fetch interval");

        // read the subsequence defined by the interval into a vector
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Couldn't read the interval");

        // convert to string
        let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");

        // Split the string by whitespace and collect into a Vec<&str>
        let alleles_list: Vec<&str> = alleles.split_whitespace().collect();

        // insert alleles into reference sequence to get variable sequences 
        // replace reference allele with variant allele
        let mut var_seqs = Vec::new();
        let ku = *k as usize;
        let ref_allele_len = alleles_list[0].len();

        for allele in &alleles_list {
            let mut var_seq = seq_string.clone();
            var_seq.replace_range(ku..ku+ref_allele_len, allele);
            var_seqs.push(var_seq);
        }

        // get k-mers
        let mut kmers_by_allele = Vec::new();
        for var_seq in &var_seqs {
            let kmer_list: Vec<String> = common::get_canonical_kmers(var_seq, ku);
            kmers_by_allele.push(kmer_list);
        }

        // remove kmers shared across alleles
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


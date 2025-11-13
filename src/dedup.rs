use std::fs::File;
use std::io::{Write, BufWriter};
use std::io::{self, prelude::*, BufReader};
use std::collections::HashSet;
use std::collections::HashMap;
use rust_htslib::bcf::Reader;
use bio::io::fasta::IndexedReader;
use rust_htslib::bcf::Read;
use bio::bio_types::genome::AbstractLocus;

use crate::common;

// look across non-variable reference regions
// remove any allele-specific k-mers also found in these regions
pub fn reference_hashset(fasta_path: &String, vcf_path: &String, k: &i64) -> HashSet<String> {
    println!("Build hashset of reference k-mers...");
    // rust-htslib provides VCF I/O.
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");

    // read indexed fasta file
    let mut faidx = IndexedReader::from_file(fasta_path).unwrap();

    // loop over variants a
    let mut start = 1;
    let mut result = Vec::new();
    for record_result in vcf_reader.records() {
        let record = record_result.expect("Failure reading record");
        let pos = record.pos() - 1;
        let chrom = record.contig();
        let end = pos - 1;
        if (end - start) <= *k {
            println!("Skipping adjacent variants at: {} {} {}", chrom, start, end);
            continue
        }
        // move the pointer in the index to the desired sequence and interval
        faidx.fetch(chrom, start.try_into().unwrap(), end.try_into().unwrap()).expect("Could not fetch interval");
        // read the subsequence defined by the interval into a vector
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Could not read interval");
        // convert to string
        let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
        // get k-mers
        let ref_kmers: Vec<String> = common::get_canonical_kmers(&seq_string, *k as usize);
        result.extend(ref_kmers);
        start = pos + 1;
    }
    // put kmers into hashset
    let ref_kmers_hashset: HashSet<String> = result.into_iter().collect();
    return ref_kmers_hashset;
}

// open index and remove k-mers found in reference
pub fn remove_ref_kmers(index: &String, output: &String, ref_hashset: HashSet<String>) -> Result<Vec<Vec<Vec<String>>>, io::Error> {

    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<String>> = kmers_list
        .iter()
        .map(|s| s.split(';').map(|x| x.to_string()).collect())
        .collect();

        data.push(kmers_by_allele);
    }

    // remove any k-mers found in reference
    for inner_vec in &mut data {
        for inner_inner_vec in inner_vec{
            inner_inner_vec.retain(|s| !ref_hashset.contains(s));
        }
    }

    // output file
    let mut buffered_file = BufWriter::new(File::create(output)?);
    // write remaining k-mers to new index
    // open index again
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut i = 0;
    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let dedup_kmers = &data[i];
        let num_kmers_per_allele = dedup_kmers.into_iter().map( |inner_vec| inner_vec.len().to_string()).collect::<Vec<String>>().join("|");
        let dedup_string = dedup_kmers.into_iter().map(|inner_vec| inner_vec.join(";")).collect::<Vec<String>>().join("|");
        let parts = vec![fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], &num_kmers_per_allele, &dedup_string];
        writeln!(buffered_file, "{}", parts.join(","))?;
        i += 1;
    }
    Ok(data)
}

// look across variants to find shared kmers
pub fn find_dup_kmers_across_var(index: &String, output: &String) -> Result<Vec<Vec<Vec<String>>>, io::Error> {
    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<String>> = kmers_list
        .iter()
        .map(|s| s.split(';').map(|x| x.to_string()).collect())
        .collect();

        data.push(kmers_by_allele);
    }

    // first pass: count how many times k-mers are found across alleles
    let mut counts: HashMap<String, usize> = HashMap::new();

    for inner_vec in &data {
        for inner_inner_vec in inner_vec {
            for s in inner_inner_vec{
                *counts.entry(s.to_string()).or_insert(0) += 1;
            }
        }
    }

    // identify k-mers found more than once
    let dup_kmers: Vec<String> = counts
        .into_iter() // Get an iterator over key-value pairs
        .filter(|(_key, value)| *value > 1) // Filter pairs where value > 1
        .map(|(key, _value)| key) // Map to get only the keys
        .collect(); // Collect the keys into a Vec

    let dup_kmers_hashset: HashSet<String> = dup_kmers.into_iter().collect();

    // second pass: remove any k-mers that were duplicated
    for inner_vec in &mut data {
        for inner_inner_vec in inner_vec{
            inner_inner_vec.retain(|s| !dup_kmers_hashset.contains(s));
        }
    }

    // rewrite index with deduped k-mers
    // output file
    let mut buffered_file = BufWriter::new(File::create(output)?);

    // loop through vector of deduped k-mers
    let mut i = 0;
    
    // open index again
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();

        let dedup_kmers = &data[i];

        let num_kmers_per_allele = dedup_kmers.into_iter().map( |inner_vec| inner_vec.len().to_string()).collect::<Vec<String>>().join("|");

        let dedup_string = dedup_kmers.into_iter().map(|inner_vec| inner_vec.join(";")).collect::<Vec<String>>().join("|");
        
        // write index to output file
        let parts = vec![fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], &num_kmers_per_allele, &dedup_string];
        writeln!(buffered_file, "{}", parts.join(","))?;

        i += 1;
    }

    Ok(data)
}

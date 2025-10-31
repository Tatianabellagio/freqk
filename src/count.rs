use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::collections::HashSet;
use std::collections::HashMap;
use fastq::{parse_path as other_parse_path, Record as OtherRecord};

use crate::common;

// get kmer length from index
pub fn k_from_index(index: &String) -> Result<i64, io::Error> {
    // Open the file
    let file = File::open(index)?;

    // Create a buffered reader
    let mut reader = BufReader::new(file);

    // Create a string to store the first line
    let mut first_line = String::new();

    // Read the first line into the string
    reader.read_line(&mut first_line)?;

    // parse kmer list
    let split_line: Vec<&str> = first_line.split(',').collect();
    let kmers = split_line[7];
    let kmers_list: Vec<&str> = kmers.split('|').collect();

    let kmers_by_allele: Vec<Vec<&str>> = kmers_list
    .iter()
    .map(|s| s.split(';').collect())
    .collect();

    // Access the first inner vector
    let first_inner_vec = &kmers_by_allele[0];

    // Access the first element of that inner vector
    let first_element = first_inner_vec[0];

    return Ok(first_element.len() as i64);
}


// loads index from file and puts all allele-specific k-mers into a single hashset
pub fn build_kmer_hashset(index: String) -> Result<HashSet<String>, io::Error> {

    println!("Loading index at {} into a hashset...", index);

    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    // HashSet to store all k-mers
    let mut kmers_hashset: HashSet<String> = HashSet::new();

    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();
        
        let kmers_by_allele: Vec<Vec<&str>> = kmers_list
        .iter()
        .map(|s| s.split(';').collect())
        .collect();
        
        let kmers_all: Vec<&str> = kmers_by_allele
        .into_iter() // Consumes the outer Vec and produces an iterator over inner Vecs
        .flatten()   // Flattens the iterator of inner Vecs into an iterator of &str
        .collect();

        // convert to owned strings
        let kmers_string: Vec<String> = kmers_all.iter().map(|s| s.to_string()).collect();

        kmers_hashset.extend(kmers_string.iter().cloned());
    }

    Ok(kmers_hashset)
}

fn merge_hashmaps(vec_of_maps: Vec<HashMap<String, usize>>) -> HashMap<String, usize> {
    let mut merged_map: HashMap<String, usize> = HashMap::new();

    for map in vec_of_maps {
        for (key, value) in map {
            *merged_map.entry(key).or_insert(0) += value;
        }
    }
    merged_map
}

pub fn count_target_kmers_in_reads(index: String, reads: &String, k: i64, nthreads: usize) -> HashMap<String, usize> {
    let merged_counts: HashMap<String, usize> = other_parse_path(Some(reads), |parser| {
        let results: Vec<HashMap<String, usize>> = parser.parallel_each(nthreads, move |record_sets| {
            let mut kmer_counts: HashMap<String, usize> = HashMap::new();
            let kmers_hashset = build_kmer_hashset(index.clone()).expect("Reason");
            let mut num_records = 0;
            println!("Looping over reads...");
            for record_set in record_sets {
                for record in record_set.iter() {
                    if num_records % 10000 == 0 {
                        println!("Reads processed: {}", num_records);
                    }
                    num_records += 1;
                    let read_kmers = common::get_canonical_kmers( std::str::from_utf8( record.seq() ).expect("REASON"), k.try_into().unwrap());
                        for read_kmer in &read_kmers{
                            if kmers_hashset.contains(read_kmer) {
                                let count = kmer_counts.entry(read_kmer.to_string()).or_insert(0);
                                *count += 1; // Increment the count
                            }
                        }
                }
            }

            // The values we return (it can be any type implementing `Send`)
            // are collected from the different threads by
            // `parser.parallel_each` and returned. See doc for a description of
            // error handling.
            kmer_counts
        }).expect("Invalid fastq file");
        println!("Merging hashmaps...");
        let merged_counts = merge_hashmaps(results);
        //println!("{:?}", kmer_counts);
        merged_counts
    }).expect("Invalid compression");
    return merged_counts;
}

// write k-mer counts to file
pub fn write_kmers(kmer_counts: HashMap<String, usize>, output: &String) -> io::Result<()> {
    println!("Writing k-mer counts to {}", output);
    // Open the file for writing
    let mut file = File::create(output)?;

    // Iterate over the HashMap and write each entry to the file
    for (key, value) in kmer_counts.iter() {
        writeln!(file, "{}\t{}", key, value)?;
    }

    println!("Counts successfully written to {}", output);

    Ok(())
}

// turn k-mer counts into totals per allele
// should reflect underlying allele frequencies
pub fn combine_counts_by_allele(index: &String, counts: HashMap<String, usize>) -> Result<Vec<String>, io::Error> {

    println!("Combining counts by allele...");

    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut result = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<&str>> = kmers_list
        .iter()
        .map(|s| s.split(';').collect())
        .collect();

        let mut counts_by_allele = Vec::new();

        for allele in kmers_by_allele{
            let mut total_allele_count: usize = 0;
            for kmer in allele{
                if let Some(kmer_count) = counts.get(kmer) {
                    total_allele_count += *kmer_count;
                }
            }
            counts_by_allele.push(total_allele_count);
        }
        
        let counts_by_allele_str = counts_by_allele.iter().map(|num| num.to_string()).collect::<Vec<String>>().join("|");
        
        result.push(counts_by_allele_str);
    }
    Ok(result)
}

pub fn write_counts_by_allele(counts_by_allele: Vec<String>, output: &String) -> io::Result<()>{

    println!("Writing counts by allele to {}", output);
    // Open the file for writing
    let mut file = File::create(output)?;

    // Iterate over the HashMap and write each entry to the file
    for element in counts_by_allele.iter() {
        writeln!(file, "{}", element)?;
    }

    println!("Counts successfully written.");

    Ok(())

}

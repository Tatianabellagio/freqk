use std::fs::File;
use std::io::{Write, BufWriter};
use std::io::{self, prelude::*, BufReader};
use std::collections::HashSet;
use std::collections::HashMap;

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

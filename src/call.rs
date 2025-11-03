use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use crate::common;

// array of allele-specific k-mer counts
struct Karray {
    counts: Vec<Vec<f32>>,
}

impl Karray {
    // write array to file
    fn write(self, output_file: &String) {
        let mut output_str = Vec::new();
        for one_site in &self.counts {
            let one_site_str = one_site.iter().map(|num| num.to_string()).collect::<Vec<String>>().join("|");
            output_str.push(one_site_str);
        }
        let _ = common::write_strings(output_str, output_file);
    }
}

// converting string array to u32 array
fn vec_vec_string_to_vec_vec_u32(string_vec_vec: Vec<Vec<String>>) -> Vec<Vec<u32>> {
    let result: Vec<Vec<u32>> = string_vec_vec
        .into_iter()
        .map(|inner_vec_strings| {
            inner_vec_strings
                .into_iter()
                .filter_map(|s| s.parse::<u32>().ok())
                .collect()
        })
        .collect();

    return result;
}

// dividing two 2d arrays
fn elementwise_division_2d(vec_a: &Vec<Vec<u32>>, vec_b: &Vec<Vec<u32>>) -> Result<Vec<Vec<f32>>, io::Error> {
    if vec_a.len() != vec_b.len() {
        panic!("Vectors of different lengths");
    }

    let result: Vec<Vec<f32>> = vec_a
        .iter()
        .zip(vec_b.iter())
        .map(|(inner_a, inner_b)| {
            if inner_a.len() != inner_b.len() {
                panic!("Inner vectors of different lengths");
            }
            inner_a
                .iter()
                .zip(inner_b.iter())
                .map(|(&a, &b)| {
                    if b == 0 {
                        panic!("Division by zero detected");
                    }
                    (a as f32)/(b as f32)
                })
                .collect()
        })
        .collect();

    Ok(result)
}

// find minimum number in Vec<u32>
fn normalized_counts_to_allele_freq(norm_counts: Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    // get total of normalized counts
    //println!("{:?}", norm_counts);

    let mut sums: Vec<f32> = Vec::new();
    for inner_vec in &norm_counts {
        sums.push(inner_vec.iter().sum());
    }

    //println!("{:?}", sums);
    let mut result = Vec::new();
    // divide normalized frequencies by sums
    for (norm_count, sum) in norm_counts.iter().zip(sums.iter()) {
        let mut allele_freqs = Vec::new();
        for count in norm_count{
            allele_freqs.push(count/sum)
        }
        result.push(allele_freqs);
    }
    return result;
}

// BRING IT ALL TOGETHER! :D
pub fn call_from_counts(index: &String, counts: &String, output: &String) -> Result<(), io::Error> {
    // parse number of allele-specific k-mers from index
    let num_uniq_kmers_per_allele = common::read_index_field(index, 6);
    // convert to u32
    let num_uniq_kmers_per_allele_u32 = vec_vec_string_to_vec_vec_u32(num_uniq_kmers_per_allele?);
    // parse counts file
    let file = File::open(counts)?;
    let reader = BufReader::new(file);
    let mut kmer_counts_per_allele = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<String> = line.split('|').map(|s| s.to_owned()).collect();
        kmer_counts_per_allele.push(fields);
    }
    // convert to u32
    let kmer_counts_per_allele_u32 = vec_vec_string_to_vec_vec_u32(kmer_counts_per_allele);
    // normalize k-mer counts by number of alleles
    let counts_per_kmer = elementwise_division_2d(&kmer_counts_per_allele_u32, &num_uniq_kmers_per_allele_u32);
    // get allele frequencies 
    let allele_freq = Karray {
        counts: normalized_counts_to_allele_freq(counts_per_kmer?),
    };

    // write frequencies to file
    allele_freq.write(output);
    Ok(())
}

use std::io;
use std::fs::File;
use std::io::Write;
use crate::common::io::BufReader;
use std::io::BufRead;

// slide window of k bases through sequence
// get reverse complements
// keep lexicographically sooner k-mer for each pair
pub fn get_canonical_kmers(sequence: &str, k: usize) -> Vec<String> {
    let mut canonical_kmers = Vec::new();

    if k == 0 || k > sequence.len() {
        return canonical_kmers; // Handle invalid k-mer length
    }

    for i in 0..=(sequence.len() - k) {
        let kmer_slice = &sequence[i..i + k];
        let reverse_complement = reverse_complement(kmer_slice);

        // Compare lexicographically to find the canonical k-mer
        if kmer_slice < reverse_complement.as_str() {
            canonical_kmers.push(kmer_slice.to_string());
        } else {
            canonical_kmers.push(reverse_complement.to_string());
        }
    }
    canonical_kmers
}

// Helper function to compute the reverse complement of a DNA sequence
fn reverse_complement(dna_sequence: &str) -> String {
    dna_sequence
        .chars()
        .rev() // Reverse the sequence
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => base, // Handle non-DNA characters or other cases
        })
        .collect()
}

// Read field in index
pub fn read_index_field(index: &String, column: usize) -> Result<Vec<Vec<String>>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut result = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let field = fields[column];
        let field_vec: Vec<String> = field.split('|').map(|s| s.to_owned()).collect();
        result.push(field_vec);
    }
    Ok(result)
}

// Write Vec<String> to a file
pub fn write_strings(strings: Vec<String>, output: &String) -> io::Result<()>{

    println!("Writing allele frequencies to: {}", output);
    // Open the file for writing
    let mut file = File::create(output)?;

    // Iterate over the HashMap and write each entry to the file
    for string in strings.iter() {
        writeln!(file, "{}", string)?;
    }

    println!("Writing successful");

    Ok(())

}


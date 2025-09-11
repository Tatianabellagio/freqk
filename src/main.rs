use clap::{Parser, Subcommand};
use rust_htslib::bcf::Reader;
use rust_htslib::bcf::Read;
use bio::io::fasta::IndexedReader;
use std::fs::File;
use std::io::{Write, BufWriter};
use bio::bio_types::genome::AbstractLocus;
use std::io::{self, prelude::*, BufReader};
use std::collections::HashSet;
use kseq::parse_path;
use std::collections::HashMap;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}


#[derive(Subcommand, Debug)]
enum Commands {
    /// Creates a new user
    Index {
        #[arg(long)]
        fasta: String,
        #[arg(long)]
        vcf: String,
        #[arg(long)]
        output: String,
        #[arg(long)]
        k: i64,
    },
    /// Deletes an existing user
    Count {
        #[arg(long)]
        index: String,
        #[arg(long)]
        reads: String,
        #[arg(long)]
        output: String,
        #[arg(long)]
        k: i64,
    },
    /// Lists all users
    Call {
        #[arg(long)]
        counts: String,
        #[arg(long)]
        output: String,
    },
}

// slide window of k bases through sequence
// get reverse complements
// keep lexicographically sooner k-mer for each pair
fn get_canonical_kmers(sequence: &str, k: usize) -> Vec<String> {
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


// insert variant into reference, get k-mers for each variant
fn insert_var(vcf_path: &String, fasta_path: &String, output_path: &String, k: &i64) -> Option<Vec<i32>> {
    // rust-htslib provides VCF I/O.
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");

    // read indexed fasta file
    let mut faidx = IndexedReader::from_file(fasta_path).unwrap();

    // output file
    let mut buffered_file = BufWriter::new(File::create(output_path).ok()?);
    //let mut file = File::create(path)?; // Creates or overwrites the file

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
        let pos = record.pos();
        let chrom = record.contig();
        faidx.fetch(chrom, (pos - k).try_into().unwrap(), (pos + k).try_into().unwrap() ).expect("Couldn't fetch interval");

        // read the subsequence defined by the interval into a vector
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Couldn't read the interval");

        // convert to string
        let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");

        // Split the string by whitespace and collect into a Vec<&str>
        let alleles_list: Vec<&str> = alleles.split_whitespace().collect();

        //let max_len = get_max_allele_length(alleles_list);

        // insert alleles into reference sequence to get variable sequences 
        let mut var_seqs = Vec::new();
        let ku = *k as usize;

        for allele in &alleles_list {
            let mut var_seq = seq_string.clone();
            var_seq.replace_range(ku..ku+alleles_list[0].len(), allele);
            var_seqs.push(var_seq);
        }

        //println!("{}", alleles_list[0]);
        //println!("{:?}", seq_string);
        //seq_string.replace_range(ku..ku+1, alleles_list[1]);
        //println!("{:?}", seq_string);
        
        // get k-mers
        let mut joined_kmers_list = Vec::new();
        for var_seq in &var_seqs {
            let kmer_list: Vec<String> = get_canonical_kmers(var_seq, ku);
            let joined_kmers = kmer_list.join(";");
            joined_kmers_list.push(joined_kmers);
        }


        // write index to output file
        let parts = vec![i.to_string(), chrom.to_string(), pos.to_string(), seq_string.clone(), alleles_list.join("|"), var_seqs.join("|"), joined_kmers_list.join("|")];
        writeln!(buffered_file, "{}", parts.join(",")).ok()?;

        // 0-based position and the list of alleles
        //println!("{}, Chrom: {}, Locus: {}, Alleles: {:?}, Var Seqs: {:?}, Seq: {:?}, kmers: {:?}", i, chrom, pos, alleles_list, var_seqs, seq_string, joined_kmers_list);
    }
    None
}

















// read index by line
fn build_kmer_hashset(index: &String) -> Result<HashSet<String>, io::Error> {

    println!("Loading index at {} into a hashset...", index);

    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    // HashSet to store all k-mers
    let mut kmers_hashset: HashSet<String> = HashSet::new();

    for line_result in reader.lines() {
        //println!("{}", &line?);
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        //let parts: Vec<&str> = line?.split(',').collect();
        let kmers = fields[6];
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

        // Extend the HashSet with elements from the vector
        kmers_hashset.extend(kmers_string.iter().cloned());
        
        //let mut kmers: Vec<&str> = kmers.split('|').collect();
        //let kmers_list: Vec<&str> = kmers.split(';').collect();
        //let kmers_list: Vec<&str> = kmers_list.split('|').collect();
        //println!("{:?}", kmers_hashset);
    }

    Ok(kmers_hashset)
}


// input hashset of kmers, path to reads
// loop over reads and get canonical kmers
// only count k-mer if its in hashset
// save output to k-mer counts file
fn count_target_kmers_in_reads(kmers_hashset: HashSet<String>, reads: &String, k: i64) -> HashMap<String, usize> {
    println!("Counting k-mers in reads...");
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
    let mut records = parse_path(reads).unwrap();
	// let mut records = parse_reader(File::open(path).unwrap()).unwrap();
	while let Some(record) = records.iter_record().unwrap() {
		//println!("head:{} des:{} seq:{} qual:{} len:{}", 
		//	record.head(), record.des(), record.seq(), 
		//	record.qual(), record.len());
        let read_kmers = get_canonical_kmers(record.seq(), k.try_into().unwrap());
        for read_kmer in &read_kmers{
            if kmers_hashset.contains(read_kmer) {
                //println!("kmer is in hashset!");
                let count = kmer_counts.entry(read_kmer.to_string()).or_insert(0);
                *count += 1; // Increment the count
            }
        }
	}
    kmer_counts
}




// write k-mer counts to file
fn write_kmers(kmer_counts: HashMap<String, usize>, output: &String) -> io::Result<()> {
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


fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Index { fasta, vcf, output, k } => {
            println!("fasta: {}, vcf: {}, k: {}", fasta, vcf, k);
            insert_var(vcf, fasta, output, k);
        }
        Commands::Count { index, reads, k, output } => {
            println!("Counting k-mers: {}, {}, {}, {}", index, reads, k, output);
            let kmer_hashset = build_kmer_hashset(index);
            let kmer_counts = count_target_kmers_in_reads(kmer_hashset.expect("Error creating target kmer hashset"), reads, *k);
            write_kmers(kmer_counts, output);
        }
        Commands::Call { counts, output } => {
            println!("Converting counts to allele frequencies... {}, {}", counts, output);
        }
    }
}

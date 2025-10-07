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
use seq_io::fastq::{Reader as OtherReader, Record};
use seq_io::parallel::parallel_fastq;
use fastq::{parse_path as other_parse_path, Record as OtherRecord};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}


#[derive(Subcommand, Debug)]
enum Commands {
    /// Get k-mers specific to each allele of each variant
    Index {
        #[arg(short,long, help = "fasta file of reference genome")]
        fasta: String,
        #[arg(short,long, help = "vcf file of variations between reference and other genomes")]
        vcf: String,
        #[arg(short,long, help = "name of the index file to be made")]
        output: String,
        #[arg(short,long, help = "kmer length")]
        kmer: i64,
    },
    /// Deduplicate index: remove k-mers shared across variants
    Dedup {
        #[arg(short,long, help = "name of index file")]
        index: String,
        #[arg(short,long, help = "name of deduplicated index")]
        output: String,
    },
    /// Count k-mers by allele
    Count {
        #[arg(short,long, help = "name of index files")]
        index: String,
        #[arg(short,long, help = "comma-separated list of reads in fastq format (can be gz or not)", value_delimiter = ',', value_name = "FILE")]
        reads: String,
        #[arg(short,long, help = "Number of threads to use (a copy of index will be loaded onto each thread)")]
        nthreads: usize,
        #[arg(short,long, help = "name of output for counts by allele")]
        freq_output: String,
        #[arg(short,long, help = "name of output for raw kmer counts")]
        count_output: String,
    },
    /// Convert k-mer counts to an allele frequency
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
        .into_iter() // Get an iterator over key-value pairs
        .filter(|(_key, value)| *value > 1) // Filter pairs where value > 2
        .map(|(key, _value)| key) // Map to get only the keys
        .collect(); // Collect the keys into a Vec

    let dup_kmers_hashset: HashSet<String> = dup_kmers.into_iter().collect();

    // second pass: remove any k-mers that were duplicated
    for inner_vec in &mut data {
        inner_vec.retain(|s| !dup_kmers_hashset.contains(s));
    }

    data
}


// look across variants to find shared kmers
fn find_dup_kmers_across_var(index: &String, output: &String) -> Result<Vec<Vec<Vec<String>>>, io::Error> {
    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();

    for line_result in reader.lines() {
        //println!("{}", &line?);
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

    //Ok(data)

    // rewrite index with deduped k-mers
    // output file
    let mut buffered_file = BufWriter::new(File::create(output)?);

    // loop through vector of deduped k-mers
    let mut i = 0;
    
    // open index again
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    for line_result in reader.lines() {
        //println!("{}", &line?);
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

        //let max_len = get_max_allele_length(alleles_list);

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
            let kmer_list: Vec<String> = get_canonical_kmers(var_seq, ku);
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


// read index by line
fn build_kmer_hashset(index: String) -> Result<HashSet<String>, io::Error> {

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

        // Extend the HashSet with elements from the vector
        kmers_hashset.extend(kmers_string.iter().cloned());
    }

    Ok(kmers_hashset)
}


// input hashset of kmers, path to reads
// loop over reads and get canonical kmers
// only count k-mer if its in hashset
// save output to k-mer counts file
fn count_target_kmers_in_reads(kmers_hashset: HashSet<String>, reads: &Vec<String>, k: i64) -> HashMap<String, usize> {
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
    for read_set in reads{
    println!("Counting k-mers in {}...", read_set);
    let mut records = parse_path(read_set).unwrap();
    let mut i = 0;
	// let mut records = parse_reader(File::open(path).unwrap()).unwrap();
	while let Some(record) = records.iter_record().unwrap() {
		//println!("head:{} des:{} seq:{} qual:{} len:{}", 
		//	record.head(), record.des(), record.seq(), 
		//	record.qual(), record.len());
        if i % 10_000 == 0 {
            println!("Reads processed: {}", i);
        }
        let read_kmers = get_canonical_kmers(record.seq(), k.try_into().unwrap());
        for read_kmer in &read_kmers{
            if kmers_hashset.contains(read_kmer) {
                //println!("kmer is in hashset!");
                let count = kmer_counts.entry(read_kmer.to_string()).or_insert(0);
                *count += 1; // Increment the count
            }
        }
        i += 1;
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


fn write_counts_by_allele(counts_by_allele: Vec<String>, output: &String) -> io::Result<()>{

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

// turn k-mer counts into totals per allele
// should reflect underlying allele frequencies
fn combine_counts_by_allele(index: &String, counts: HashMap<String, usize>) -> Result<Vec<String>, io::Error> {

    println!("Combining counts by allele...");

    // open index
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut result = Vec::new();

    for line_result in reader.lines() {
        //println!("{}", &line?);
        let line = line_result?; // Handle potential errors reading the line
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<&str>> = kmers_list
        .iter()
        .map(|s| s.split(';').collect())
        .collect();

        //println!("{:?}", kmers_by_allele);

        let mut counts_by_allele = Vec::new();

        for allele in kmers_by_allele{
            let mut total_allele_count: usize = 0;
            for kmer in allele{
                if let Some(kmer_count) = counts.get(kmer) {
                    total_allele_count += *kmer_count;
                }
            }
            //println!("Next allele...");
            counts_by_allele.push(total_allele_count);
        }
        
        let counts_by_allele_str = counts_by_allele.iter().map(|num| num.to_string()).collect::<Vec<String>>().join("|");
        
        result.push(counts_by_allele_str);
    }
    Ok(result)
}


// get kmer length from index
fn k_from_index(index: &String) -> Result<i64, io::Error> {
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


fn count_target_kmers_in_reads_2(kmers_hashset: HashSet<String>, reads: &String, k: i64) -> HashMap<String, usize>{
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
    let reader = OtherReader::from_path(reads).unwrap();
    let mut writer = BufWriter::new(File::create("filtered.fastq").unwrap());

    let test_counts = parallel_fastq(reader, 4, 2,
        |record, found| { // runs in worker
            *found = record.seq().windows(3).position(|s| s == b"AAA").is_some();
        },
        |record, found| { // runs in main thread
            if *found {
                record.write(&mut writer).unwrap();
            }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed,
        // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
            None::<()>
        }).unwrap();

    return kmer_counts;
}

// use fastq-rs to loop over reads with multiple threads

fn merge_hashmaps(vec_of_maps: Vec<HashMap<String, usize>>) -> HashMap<String, usize> {
    let mut merged_map: HashMap<String, usize> = HashMap::new();

    for map in vec_of_maps {
        for (key, value) in map {
            *merged_map.entry(key).or_insert(0) += value;
        }
    }
    merged_map
}

fn count_target_kmers_in_reads_3(index: String, reads: &String, k: i64, nthreads: usize) -> HashMap<String, usize> {
    // Treat "-" as stdin
    //let path = match <std::string::String as AsRef<T>>::as_ref(reads).map(String::as_ref) {
    //    None | Some("-") => { None },
    //    Some(name) => Some(name)
    //};

    //let mut kmer_counts: HashMap<String, usize> = HashMap::new();

    let merged_counts: HashMap<String, usize> = other_parse_path(Some(reads), |parser| {
        let results: Vec<HashMap<String, usize>> = parser.parallel_each(nthreads, move |record_sets| {
            let mut thread_total = 0;
            let mut kmer_counts: HashMap<String, usize> = HashMap::new();
            let kmers_hashset = build_kmer_hashset(index.clone()).expect("Reason");
            //let kmers_hashset: HashSet<String> = HashSet::new();
            //let k = 31;
            println!("Looping over reads...");
            for record_set in record_sets {
                for record in record_set.iter() {
                    thread_total += 1;
                    let read_kmers = get_canonical_kmers( std::str::from_utf8( record.seq() ).expect("REASON"), k.try_into().unwrap());
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

// bring it all together!
fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Index { fasta, vcf, output, kmer } => {
            println!("fasta: {}, vcf: {}, k: {}", fasta, vcf, kmer);
            insert_var(vcf, fasta, output, kmer);
        }
        Commands::Count { index, reads, nthreads, freq_output, count_output } => {
            println!("Counting k-mers: {}, {:?}, {}, {}, {}", index, reads, nthreads, freq_output, count_output);
            let k = k_from_index(index);
            println!("k is: {:?}", k);
            //let kmer_hashset = build_kmer_hashset(index);
            let kmer_counts =  count_target_kmers_in_reads_3(index.clone(), reads, k.expect("Cannot parse kmer length from index."), *nthreads);
            //let kmer_counts = count_target_kmers_in_reads(kmer_hashset.expect("Error creating target kmer hashset"), reads, k.expect("Cannot parse kmer length from index."));
            //let _ = write_kmers(kmer_counts.clone(), count_output);
            let counts_by_allele = combine_counts_by_allele(index, kmer_counts);
            //println!("{:?}", counts_by_allele);
            let _ = write_counts_by_allele(counts_by_allele.expect("Error writing counts by allele"), freq_output);
        }
        Commands::Dedup { index, output } => {
            println!("Deduplicating index: {} {}", index, output);
            let _ = find_dup_kmers_across_var(index, output);
            //println!("{:?}", deduped_kmers_by_allele);
        }
        Commands::Call { counts, output } => {
            println!("Converting counts to allele frequencies... {}, {}", counts, output);
        }
    }
}

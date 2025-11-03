use clap::{Parser, Subcommand};

mod common;
mod index;
mod dedup;
mod count;
mod call;

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
        #[arg(short,long, help = "name of the index file to be output")]
        output: String,
        #[arg(short,long, help = "kmer length for building the index")]
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
        #[arg(short,long, help = "name of index file")]
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
    /// Convert counts by allele into allele frequencies
    Call {
        #[arg(short,long, help = "name of index file")]
        index: String,
        #[arg(short,long, help = "name of file for counts by allele")]
        counts: String,
        #[arg(short,long, help = "name of output file with allele frequencies")]
        output: String,
    },
}

// bring it all together!
fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Index { fasta, vcf, output, kmer } => {
            println!("fasta: {}, vcf: {}, k: {}", fasta, vcf, kmer);
            index::index_workflow(vcf, fasta, output, kmer);
        }
        Commands::Count { index, reads, nthreads, freq_output, count_output } => {
            println!("Counting k-mers: {}, {:?}, {}, {}, {}", index, reads, nthreads, freq_output, count_output);
            count::count_workflow(index, reads, *nthreads, freq_output, count_output);
        }
        Commands::Dedup { index, output } => {
            println!("Deduplicating index: {} {}", index, output);
            let _ = dedup::find_dup_kmers_across_var(index, output);
        }
        Commands::Call { index, counts, output} => {
            println!("Converting counts to allele frequencies: {} {} {}", index, counts, output);
            let _ = call::call_from_counts(index, counts, output);
        }
    }
}

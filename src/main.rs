use clap::{Parser, Subcommand};
use rust_htslib::bcf::Reader;
use rust_htslib::bcf::Read;

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
        k: u32,
    },
    /// Deletes an existing user
    Count {
        #[arg(long)]
        id: u32,
    },
    /// Lists all users
    Call {
        #[arg(long)]
        id: u32,
    },
}


fn insert_var(vcf_path: &String) -> () {
    // rust-htslib provides VCF I/O.
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");

    // iterate through each row of the vcf body.
    for (i, record_result) in vcf_reader.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut s = String::new();
        for allele in record.alleles() {
            for c in allele {
                s.push(char::from(*c))
             }
            s.push(' ')
        }
    // 0-based position and the list of alleles
    println!("{}, Locus: {}, Alleles: {}", i, record.pos(), s);
    }
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Index { fasta, vcf, k } => {
            println!("fasta: {}, vcf: {}, k: {}", fasta, vcf, k);
            insert_var(vcf);
        }
        Commands::Count { id } => {
            println!("Deleting user with ID: {}", id);
            // Logic to delete a user
        }
        Commands::Call { id } => {
            println!("Listing all users... {}", id);
            // Logic to list users
        }
    }
}

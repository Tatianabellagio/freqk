// functions to check that input is properly formated
    // check that k-mers are lexicographically sorted
    pub fn check_sort(seqs: &Vec<String>) -> bool {
        // Limit to the first 1000 elements (or fewer)
        let limit = seqs.len().min(1000);
        let seqs_sub = &seqs[..limit];

        // Make a new Vec and sort it
        let mut sorted_seqs = seqs_sub.to_vec();
        sorted_seqs.sort();

        // Compare original slice with sorted one
        let result = seqs_sub == &sorted_seqs;

        println!("Input sorted: {}", result);
        result
    }

    // check that only ATGC are in alphabet
    pub fn check_letters(seqs: &Vec<String>) -> bool {
        // Limit to the first 1000 elements (or fewer)
        let limit = seqs.len().min(1000);
        let seqs_sub = &seqs[..limit];

        let result = seqs_sub.iter().all(|seq| seq.chars().all(|c| matches!(c, 'A' | 'T' | 'G' | 'C')));
        println!("Only ATGC: {}", result);
        result
    }

    // combine all checks together
    pub fn all_checks(seqs: &Vec<String>){
        println!("Checking input format...");
        if check_sort(seqs) == false {
            panic!(":(");
        }

        if check_letters(seqs) == false {
            panic!(":(")
        }
    }
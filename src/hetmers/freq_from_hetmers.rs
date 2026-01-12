use std::f64;
use statrs::function::gamma::gamma;
use average::Mean;
//use statrs::statistics::MeanN;

// empirical allele frequencies
pub fn counts_to_frequencies(count_pairs: &Vec<String>) -> Vec<String> {
    println!("Calculating frequencies...");
    let frequencies: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let min_num = num1.min(num2);
                    let max_num = num1.max(num2);
                    let sum = min_num + max_num;
                    if sum != 0.0 {
                        return Some(min_num / sum);
                    }
                }
            }
        None
        })
        .collect();

    // reformatting frequency list
    let freq_strings: Vec<String> = frequencies.into_iter().map(|s| s.to_string()).collect();
    return freq_strings;
}

// Tag hetmers with really high total coverage (potentially due to paralogous sequences)
pub fn high_cov_hetmers(count_pairs: &Vec<String>, sigma: f64, n: i32, cov: f64) -> Vec<String> {
    println!("Checking for questionable hetmers...");
    let potential_filter: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let sum = num1 + num2;
                    let stderr = ((n as f64)*(cov as f64)).sqrt();
                    if sum > (sigma*stderr) as f64 {
                        return Some(1);
                    } else{
                        return Some(0);
                    }
                }
            }
        None
        })
        .collect();

    // reformatting frequency list
    let potent_strings: Vec<String> = potential_filter.into_iter().map(|s| s.to_string()).collect();
    return potent_strings;

}

// Compute truncation constant
pub fn truncation_constant(c: usize, lambda: f64) -> f64 {
    let sum: f64 = (0..c)
        .map(|x| {
            let numerator = (-lambda).exp() * lambda.powi(x as i32);
            let denominator = gamma((x + 1) as f64); // factorial(x)
            numerator / denominator
        })
        .sum();

    1.0 - sum
}

// Compute posterior where the minor k-mer count is variable
pub fn posterior_min_kmer_count(x: f64, z: f64, n: i32, cov: f64, c: usize, alpha: f64, beta: f64) -> usize {
    let mut likelihood_times_prior = Vec::new();
    let mut total_probability = 0.0;
    let max_minor_count = n/2; // minor allele can't have frequency above 1/2 by definition

    for i in 1..max_minor_count {
        let p = i as f64 / n as f64;
        let lambda_x = (i as f64) * (cov as f64);
        let lambda_y = ((n - i) as f64) * (cov as f64);

        let tx = truncation_constant(c, lambda_x);
        let ty = truncation_constant(c, lambda_y);

        let likelihood = p.powf(x as f64 + alpha - 1.0) * (1.0 - p).powf((z - x) as f64 + beta - 1.0);
        let likelihood_truncated = likelihood / (tx * ty);
        likelihood_times_prior.push(likelihood_truncated);
        total_probability += likelihood;
    }

    let posterior: Vec<f64> = likelihood_times_prior
        .iter()
        .map(|val| val / total_probability)
        .collect();

    let max_index = posterior
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    max_index + c
}

// calculate posterior distribution for allele count
pub fn counts_to_bayes_state(count_pairs: &Vec<String>, n: i32, cov:f64, c: usize, alpha: f64, beta: f64) -> Vec<usize>{
    println!("Calculating posterior...");
    let bayes_states: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let x = num1.min(num2);
                    //let max_num = num1.max(num2);
                    //let z = x + max_num;
                    let z = num1 + num2;
                    let post = posterior_min_kmer_count(x,z,n,cov,c,alpha,beta);
                    return Some(post);
                    //return highest_prob_index(post);
                } else { return Some(0);}
            } else { return Some(0); }
        })
        .collect();
    return bayes_states;
}

pub fn freq_to_pi(freq_vec: Vec<String>) -> f64 {
    let pi_result: Mean = freq_vec.into_iter().map(|x| {
        let p = x.parse::<f64>();
        let q = 1.0 - p.clone().expect("Reason");
        1.0 - ( (p.expect("Reason")).powi(2) + (q).powi(2) )
        }).collect();
    //return average::Mean(&pi_result);
    return pi_result.mean();
}

// test functions
#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn a_few_freqs(){
        // Vec<String>, output: String
        let count_pairs = vec!["5,120".to_string(), "1,9".to_string(), "20,140".to_string(), "22,22".to_string()];
        let result = counts_to_frequencies(&count_pairs);
        let expected = vec!["0.04".to_string(), "0.1".to_string(), "0.125".to_string(), "0.5".to_string()];
        assert_eq!(result, expected);
    }
    
    // helper function to test equality of floating point numbers
    fn round_to_decimals(x: f64, decimals: u32) -> f64 {
        let factor = 10f64.powi(decimals as i32);
        (x * factor).round() / factor
    }

    #[test]
    fn small_truncation(){
        let result = round_to_decimals(truncation_constant(5, 10.0), 7);
        let expected = 0.9707473;
        assert_eq!(result, expected);
    }

    #[test]
    fn heterozygosity(){
        let input = vec!["0.04".to_string(), "0.1".to_string()];
        let expected = 0.12839999999999996;
        let result = freq_to_pi(input);
        assert_eq!(result, expected);
    }
}

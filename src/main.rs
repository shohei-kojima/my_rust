// Author: Shohei Kojima @ RIKEN

use rust_codes::io::*;

fn main() {
    let path = "/home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa".to_string();
    let mut fa = FastaRecords::from_file(&path);
    fa.to_upper();
    println!("{}", fa.fai.chrs[0]);
}

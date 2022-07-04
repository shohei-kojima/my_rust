// Author: Shohei Kojima @ RIKEN

//! IO for fasta, fai, fq

use std::collections::{HashMap, HashSet};
use std::io::BufReader;
use std::io::prelude::*;
use std::fs::File;
use std::path::PathBuf;

use crate::nucl::*;


/// Stores one line of fasta index
#[derive(Clone, Debug, PartialEq)]
pub struct FaiRecord {
    pub name:       String,
    pub length:     u64,
    pub offset:     u64,
    pub linebases:  u64,
    pub linewidth:  u64,
    pub qualoffset: u64,
    pub is_fq:      bool,
}

/// Stores fasta records
#[derive(Clone, Debug, PartialEq)]
pub struct FaiRecords {
    pub chrs: Vec<String>,
    pub map:  HashMap<String, FaiRecord>,
}

impl FaiRecord {
    /// For fai of fasta
    pub fn new(name: String, length: u64, offset: u64, linebases: u64, linewidth: u64) -> Self {
        FaiRecord{name: name, length: length, offset: offset, 
                  linebases: linebases, linewidth: linewidth, qualoffset: 0, is_fq: false}
    }
    
    /// For fai of fastq
    pub fn new_fq(name: String, length: u64, offset: u64, linebases: u64, linewidth: u64, qualoffset: u64) -> Self {
        FaiRecord{name: name, length: length, offset: offset, 
                  linebases: linebases, linewidth: linewidth, qualoffset: qualoffset, is_fq: true}
    }
}

impl FaiRecords {
    pub fn new() -> Self {
        FaiRecords{chrs: Vec::new(), map: HashMap::new()}
    }
    
    pub fn push(&mut self, s: FaiRecord) {
        if ! self.map.contains_key(&(s.name)) {
            self.chrs.push(s.name.clone());
            self.map.insert(s.name.clone(), s);
        } else {
            panic!("Adding the same chromosome ({}) in FaiRecords.", s.name);
        }
    }
    
    /// read a .fai file and returns FaiRecords
    pub fn from_file(fname: &str) -> Self {
        let mut fai = FaiRecords{chrs: Vec::new(), map: HashMap::new()};
        let fpath = PathBuf::from(fname);
        let f = File::open(fpath).expect("fai file not found.");
        for line in BufReader::new(f).lines() {
            match line {
                Ok(line) => {
                    let ls = line.trim().split("\t").collect::<Vec<&str>>();
                    if ls.len() == 5 {
                        let r = FaiRecord::new(ls[0].to_string(), 
                                               ls[1].parse::<u64>().unwrap(),
                                               ls[2].parse::<u64>().unwrap(), 
                                               ls[3].parse::<u64>().unwrap(), 
                                               ls[4].parse::<u64>().unwrap());
                        fai.push(r);
                    } else if ls.len() == 6 {
                        let r = FaiRecord::new_fq(ls[0].to_string(), 
                                                  ls[1].parse::<u64>().unwrap(),
                                                  ls[2].parse::<u64>().unwrap(), 
                                                  ls[3].parse::<u64>().unwrap(), 
                                                  ls[4].parse::<u64>().unwrap(), 
                                                  ls[5].parse::<u64>().unwrap());
                        fai.push(r);
                    } else {
                        panic!("Column number in the fai file is not 5 or 6. Please check the file format again.");
                    }
                },
                Err(e) => { panic!("{}", e); }
            }
        }
        fai
    }
    
    /// check whether FaiRecords contains a given chr
    pub fn contains(self, chr: &str) -> bool {
        self.map.contains_key(chr)
    }
}



#[derive(Clone,Debug,PartialEq)]
pub enum SeqCase {
    Upper,
    Lower,
    Mixed,
    Unknown,
}


/// Fasta reader; store all sequence in memory
#[derive(Clone, Debug, PartialEq)]
pub struct FastaRecords {
    pub map: HashMap<String, String>,
    pub attr: HashMap<String, String>,  // header attributes
    pub fai: FaiRecords,
    pub case: SeqCase,
}

impl FastaRecords {
    pub fn new() -> Self {
        FastaRecords{map: HashMap::new(), attr: HashMap::new(), fai: FaiRecords::new(), 
                     case: SeqCase::Unknown}
    }
    
    pub fn push(&mut self, chr: String, seq: String) {
        if ! self.map.contains_key(&chr) {
            self.map.insert(chr, seq);
        } else {
            panic!("The same chr already exists in FastaRecords.");
        }
    }
    
    pub fn to_upper(&mut self) {
        for (_, value) in &mut(self.map) {
            crate::nucl::to_upper_seq(value); 
        }
        self.case = SeqCase::Upper;
    }
    
    pub fn to_lower(&mut self) {
        for (_, value) in &mut(self.map) {
            crate::nucl::to_lower_seq(value); 
        }
        self.case = SeqCase::Lower;
    }
    
    pub fn to_rev_comp(&mut self) {
        for (_, value) in &mut(self.map) {
            crate::nucl::to_rev_comp_seq(value); 
        }
    }
    
    /// Takes one chr name and returns &str of the chr sequence
    pub fn get_seq(&self, chr: &str) -> Option<&str> {
        match self.map.get(chr) {
            Some(seq) => Some(&seq),
            None => None
        }
    }
    
    /// This reads fasta file with fai file
    pub fn from_file(fpath: &str) -> FastaRecords {
        let fapath = PathBuf::from(fpath);
        if ! fapath.exists() { panic!("fasta file does not exist."); }
        // make structs for fasta and fai
        let fai = FaiRecords::from_file(&(fpath.to_string() + ".fai"));
        let mut fasta = FastaRecords{map: HashMap::new(), attr: HashMap::new(), 
                                     fai: fai, case: SeqCase::Unknown};
        for chr in &(fasta.fai.chrs) {
            let mut seq = String::new();
            seq.reserve(fasta.fai.map[chr].length as usize);
            fasta.map.insert(chr.clone(), seq);
        }
        // read fasta
        let f = File::open(fpath).expect("fai file not found.");
        let mut chr: String;
        let mut attr: String;
        let mut curr_chr = "".to_string();
        for line in BufReader::new(f).lines() {
            match line {
                Ok(mut line) => {
                    if line.starts_with('>') {
                        line = line.trim()[1..].to_string();  // delete '>'
                        let mut offset = 0;
                        for (n, c) in line.char_indices() {
                            if c == ' ' || c == '\t' {
                                offset = n;
                                break;    
                            }
                        }
                        if offset > 0 {  // if attribution is present
                            chr =  (&line[.. offset]).to_string();
                            attr = (&line[offset + 1 ..]).to_string();
                        } else {
                            chr = line;
                            attr = "".to_string();
                        }
                        fasta.attr.insert(chr.clone(), attr);  // save attributions
                        curr_chr = chr;
                        if ! fasta.map.contains_key(&curr_chr) {
                            panic!("{} does not found in fai", curr_chr);
                        }
                    } else {
                        match fasta.map.get_mut(&curr_chr) {
                            Some(seq) => { seq.push_str(line.trim()); },
                            None => { panic!("{} does not found in fai", curr_chr); }
                        }
                    }
                },
                Err(e) => { panic!("{}", e); }
            }
        }
        fasta
    }
    
    /// This reads fasta file WITHOUT fai file
    pub fn from_file_wo_fai(fpath: &str) -> FastaRecords {
        let fapath = PathBuf::from(fpath);
        if ! fapath.exists() { panic!("fasta file does not exist."); }
        // make structs for fasta and fai
        let mut fasta = FastaRecords{map: HashMap::new(), attr: HashMap::new(), 
                                     fai: FaiRecords::new(), case: SeqCase::Unknown};
        // read fasta
        let f = File::open(fpath).expect("fai file not found.");
        let mut chr: String;
        let mut attr: String;
        let mut curr_chr = "".to_string();
        for line in BufReader::new(f).lines() {
            match line {
                Ok(mut line) => {
                    if line.starts_with('>') {
                        line = line.trim()[1..].to_string();  // delete '>'
                        let mut offset = 0;
                        for (n, c) in line.char_indices() {
                            if c == ' ' || c == '\t' {
                                offset = n;
                                break;    
                            }
                        }
                        if offset > 0 {  // if attribution is present
                            chr =  (&line[.. offset]).to_string();
                            attr = (&line[offset + 1 ..]).to_string();
                        } else {
                            chr = line;
                            attr = "".to_string();
                        }
                        let mut seq = String::new();
                        fasta.map.insert(chr.clone(), seq);  // add seq String
                        fasta.attr.insert(chr.clone(), attr);  // save attributions
                        fasta.fai.chrs.push(chr.clone());  // add chr in fai
                        curr_chr = chr;
                    } else {
                        match fasta.map.get_mut(&curr_chr) {
                            Some(seq) => { seq.push_str(line.trim()); },
                            None => { panic!("{} does not found in fai", curr_chr); }
                        }
                    }
                },
                Err(e) => { panic!("{}", e); }
            }
        }
        // fill fai information
        for chr in &(fasta.fai.chrs) {
            let mut r = FaiRecord::new(chr.clone().to_string(), fasta.map[chr].len() as u64, 0, 0, 0);
            fasta.fai.map.insert(chr.clone().to_string(), r);
        }
        fasta
    }
}




#[cfg(test)]
mod tests {
    use crate::io::*;
    
    #[test]
    fn test_fai() {
        let path = "./test_files/test.fa.fai".to_string();
        let fai = FaiRecords::from_file(&path);
        let mut n = 0;
        for chr in &(fai.chrs) {
            println!("{:?}", fai.map[chr]);
            n += 1;
            if n == 5 { break; }
        }
    }
    
    #[test]
    fn test_load_fasta() {
        let path = "./test_files/test.fa".to_string();
        let fa = FastaRecords::from_file(&path);
        println!("{:?}", fa);
        println!("seq of chr1: {}", fa.get_seq("chr1").unwrap());
        
        let path = "./test_files/test.fa".to_string();
        let mut fa = FastaRecords::from_file_wo_fai(&path);
        fa.to_upper();
        println!("{:?}", fa);
        println!("seq of chr1: {}", fa.get_seq("chr1").unwrap());
    }
} // mod tests

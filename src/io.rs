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



/// Fasta reader; store all sequence in memory; read with fai
pub struct FastaRecords {
    pub map: HashMap<String, String>,
    pub fai: FaiRecords,
    pub is_upper: bool,
}

impl FastaRecords {
    pub fn new() -> Self {
        FastaRecords{map: HashMap::new(), fai: FaiRecords::new(), is_upper: false}
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
    }
    
    pub fn to_lower(&mut self) {
        for (_, value) in &mut(self.map) {
            crate::nucl::to_lower_seq(value); 
        }
    }
    
    pub fn to_rev_comp(&mut self) {
        for (_, value) in &mut(self.map) {
            crate::nucl::to_rev_comp_seq(value); 
        }
    }
}




#[cfg(test)]
mod tests {
    use crate::io::*;
    
    #[test]
    fn test_fai() {
        let path = "./test_files/fasta_index.fai".to_string();
        let fai = FaiRecords::from_file(&path);
        let mut n = 0;
        for chr in &(fai.chrs) {
            println!("{:?}", fai.map[chr]);
            n += 1;
            if n == 5 { break; }
        }
    }
} // mod tests

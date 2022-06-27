//! Build interval tree and perform insersect. 

use bio::data_structures::interval_tree::{IntervalTree};
use std::collections::{HashMap, HashSet};

static INTERSECT_KEEP_DUPLICATE: bool = true;
static INTERSECT_DROP_DUPLICATE: bool = false;


/// Stores one bed record.
#[derive(Clone, Debug, PartialEq)]
pub struct BedRecord<T> {
    pub chr:    String,
    pub start:  u64,
    pub end:    u64,
    pub attr:   T,
}

impl<T> BedRecord<T> {
    pub fn new<I: Into<u64>>(chr: String, start: I, end: I, attr: T) -> Self {
        BedRecord{chr: chr, start: start.into(), end: end.into(), attr: attr}
    }
}

/// Stores multiple BedRecord objects.
#[derive(Clone, Debug, PartialEq)]
pub struct BedRecords<T> {
    pub map: HashMap<String, Vec<BedRecord<T>>>,  // (name of chr, vec)
    pub chrs: Vec<String>,
    pub is_sorted: HashMap<String, bool>,
    pub sorted_chrs: Vec<String>,
    pub is_chr_sorted: bool,
}

impl<T> BedRecords<T> {
    pub fn new() -> Self {
        BedRecords{map: HashMap::<String, Vec<BedRecord<T>>>::new(), 
                   chrs: Vec::<String>::new(), 
                   sorted_chrs: Vec::<String>::new(),
                   is_sorted: HashMap::<String, bool>::new(),
                   is_chr_sorted: false}
    }
    
    /// Add one bed record
    pub fn push(&mut self, r: BedRecord<T>) {
        let chr = r.chr.clone();
        match self.map.get_mut(&(r.chr)) {
            Some(v) => v.push(r),
            None => {
                let mut v = Vec::<BedRecord<T>>::new();
                v.push(r);
                self.map.insert(chr.clone(), v);
                self.is_sorted.insert(chr.clone(), false);
                self.chrs.push(chr.clone());
                if self.is_chr_sorted {
                    self.is_chr_sorted = false;  // when new chr name is added, make is_chr_sorted false
                }
            }
        }
        match self.is_sorted.get_mut(&chr) {
            Some(b) => { *b = false; },  // when new record is added, make is_sorted false
            None => {}
        }
    }
    
    /// Sort a vector containing chromosome names
    pub fn sort_chrnames(&mut self) {
        if ! self.is_chr_sorted {
            self.sorted_chrs = self.chrs.to_vec();
            self.sorted_chrs.sort();
            self.is_chr_sorted = true;
        }
    }
    
    /// Sort a vector containing intervals for a given chromosome
    pub fn sort_chr(&mut self, chr: &String) {
        match self.map.get_mut(chr) {
            Some(v) => {
                v.sort_by(|a, b|
                    if a.start == b.start {
                        a.end.cmp(&(b.end))
                    } else {
                        a.start.cmp(&(b.start))
                    }
                );
                match self.is_sorted.get_mut(chr) {
                    Some(b) => { *b = true; },  // sorted
                    None => {}
                }
            },
            None => {},
        }
    }
    
    /// Sort vectors containing intervals for all chromosomes
    pub fn sort_allchr(&mut self) {
        for chr in &(self.chrs) {  
            match self.map.get_mut(chr) {
                Some(v) => {
                    v.sort_by(|a, b|
                        if a.start == b.start {
                            a.end.cmp(&(b.end))
                        } else {
                            a.start.cmp(&(b.start))
                        }
                    );
                    match self.is_sorted.get_mut(chr) {
                        Some(b) => { *b = true; },  // sorted
                        None => {}
                    }
                },
                None => {},
            }
        }
    }
}


/// Similar func as `bedtools intersect -wa`
pub fn intersect_wa<'a, T1, T2>(a: &'a mut BedRecords<T1>, 
                                b: &'a mut BedRecords<T2>, 
                                allow_duplicate: bool) 
        -> (Vec<String>, HashMap::<String, Vec<&'a BedRecord<T1>>>) {
    // make a vector
    let mut intersect_map: HashMap<String, Vec<&'a BedRecord<T1>>> = HashMap::new();
    let mut intersect_chrs: Vec<String> = Vec::new();

    // Construct tree
    let mut tree;
    let mut n: u64;  // serial numbering for elements in the tree
    let mut set: HashSet<u64>;  // keeps serial numbers that are added to the vec in HashMap
    
    if ! a.is_chr_sorted { a.sort_chrnames(); }
    if ! b.is_chr_sorted { b.sort_chrnames(); }
    for chr in &(a.sorted_chrs) {
        // if not chr in b, skip
        if ! b.map.contains_key(chr) { continue; }
        // if not b is sorted, sort
        match b.is_sorted.get(chr) {
            Some(i) => {
                if ! i { b.sort_chr(chr); }
            },
            None => {}
        }
        
        // initialize tree, serial numbering, set
        tree = IntervalTree::new();
        n = 0;
        set = HashSet::new();
        for r in &(a.map[chr]) {
            tree.insert(r.start .. r.end, (r, n));
            n += 1;
        }
        for r in &(b.map[chr]) {
            for i in tree.find(r.start .. r.end) {
                match intersect_map.get_mut(chr) {
                    Some(v) => {
                        if allow_duplicate {
                            v.push(i.data().0);
                        } else if ! set.contains(&(i.data().1)) {
                            v.push(i.data().0);
                            set.insert(i.data().1);
                        }
                    },
                    None => {
                        let mut v: Vec<&'a BedRecord<T1>> = Vec::new();
                        v.push(i.data().0);
                        if ! allow_duplicate { set.insert(i.data().1); }
                        intersect_map.insert(chr.to_string(), v);
                        intersect_chrs.push(chr.to_string());
                    }
                }
            }
        }
    }
    (intersect_chrs, intersect_map)
}


/// Similar func as `bedtools intersect -wa -wb`
pub fn intersect_wawb<'a, T1, T2>(a: &'a mut BedRecords<T1>, 
                                  b: &'a mut BedRecords<T2>) 
        -> (Vec<String>, HashMap::<String, Vec<(&'a BedRecord<T1>, &'a BedRecord<T2>)>>) {
    // make a vector
    let mut intersect_map: HashMap<String, Vec<(&'a BedRecord<T1>, &'a BedRecord<T2>)>> = HashMap::new();
    let mut intersect_chrs: Vec<String> = Vec::new();
    
    // construct tree
    let mut tree;
    
    if ! a.is_chr_sorted { a.sort_chrnames(); }
    if ! b.is_chr_sorted { b.sort_chrnames(); }
    for chr in &(a.sorted_chrs) {
        if ! b.map.contains_key(chr) { continue; }
        tree = IntervalTree::new();
        for r in &(a.map[chr]) {
            tree.insert(r.start .. r.end, r);
        }
        for r in &(b.map[chr]) {
            for i in tree.find(r.start .. r.end) {
                match intersect_map.get_mut(chr) {
                    Some(v) => v.push((i.data(), r)),
                    None => {
                        let mut v = Vec::<(&'a BedRecord<T1>, &'a BedRecord<T2>)>::new();
                        v.push((i.data(), r));
                        intersect_map.insert(chr.to_string(), v);
                        intersect_chrs.push(chr.to_string());
                    }
                }
            }
        }
    }
    (intersect_chrs, intersect_map)
}


/// Similar func as `bedtools intersect -v`
pub fn intersect_v<'a, T1, T2>(a: &'a mut BedRecords<T1>, 
                               b: &'a mut BedRecords<T2>) 
        -> (Vec<String>, HashMap::<String, Vec<&'a BedRecord<T1>>>) {
    // make a vector
    let mut intersect_map: HashMap<String, Vec<&'a BedRecord<T1>>> = HashMap::new();
    let mut intersect_chrs: Vec<String> = Vec::new();
    
    // construct tree
    let mut tree;
    let mut n: u64;
    let mut v: Vec<&'a BedRecord<T1>>;
    
    // sort
    if ! a.is_chr_sorted { a.sort_chrnames(); }
    if ! b.is_chr_sorted { b.sort_chrnames(); }
    a.sort_allchr();  // if not a is sorted, sort
    
    // intersect v
    for chr in &(a.sorted_chrs) {
        v = Vec::new();
        if ! b.map.contains_key(chr) {
            for r in &(a.map[chr]) {
                v.push(r);
            }
        } else {
            // build tree
            tree = IntervalTree::new();
            for r in &(b.map[chr]) {
                tree.insert(r.start .. r.end, r);
            }
            for r in &(a.map[chr]) {
                n=0;
                for _ in tree.find(r.start .. r.end) {
                    n += 1;
                }
                if n == 0 {
                    v.push(r);
                }
            }
        }
        if v.len() >= 1 {
            intersect_map.insert(chr.to_string(), v);
            intersect_chrs.push(chr.to_string());
        }
    }
    (intersect_chrs, intersect_map)
}


/// Similar func as `bedtools merge`
pub fn merge<T>(a: &mut BedRecords<T>) -> BedRecords<()> {
    let mut merged = BedRecords::new();
    let mut s;
    let mut e;
    
    a.sort_allchr();  // if not a is sorted, sort    
    for chr in &(a.chrs) {
        if a.map[chr].len() == 0 {
            continue;
        }
        
        s = a.map[chr][0].start.clone();
        e = a.map[chr][0].end.clone();
        for r in &(a.map[chr]) {
            if r.start > e {
                merged.push(BedRecord::new(chr.clone(), s, e, ()));
                s = r.start.clone();
                e = r.end.clone();
            } else {
                e = r.end.clone();
            }
        }
        merged.push(BedRecord::new(chr.clone(), s, e, ()));
    }
    if merged.chrs.len() > 0 {
        merged.sort_chrnames();
    }
    merged
}


#[inline]
fn get_larger(x: &u64, y: &u64) -> u64 {
    if *x >= *y { (*x).clone() }
    else        { (*y).clone() }
}

#[inline]
fn get_smaller(x: &u64, y: &u64) -> u64 {
    if *x < *y { (*x).clone() }
    else       { (*y).clone() }
}

#[derive(Clone,Debug,PartialEq,Eq)]
struct SimpleRange {
    s: u64,  // start
    e: u64,  // end
}


/// Similar func as simple `bedtools intersect -a input1.bed -b input2.bed`
pub fn intersect_a<'a, T1: Clone, T2>(a: &'a mut BedRecords<T1>, 
                                      b: &'a mut BedRecords<T2>) 
        -> BedRecords<T1> {
    // make a vector
    let mut new_bed = BedRecords::new();
    
    // Construct tree
    let mut tree;
    let mut tmp_v: Vec<SimpleRange> = Vec::new(); 
    let mut n: u64;
    
    if ! a.is_chr_sorted { a.sort_chrnames(); }
    if ! b.is_chr_sorted { b.sort_chrnames(); }
    for chr in &(b.sorted_chrs) {
        // if not chr in a, skip
        if ! a.map.contains_key(chr) { continue; }
        
        // initialize tree, serial numbering, set
        tree = IntervalTree::new();
        for r in &(b.map[chr]) {
            tree.insert(r.start .. r.end, ());
        }
        for r in &(a.map[chr]) {
            n = 0;
            for i in tree.find(r.start .. r.end) {
                tmp_v.push(SimpleRange{ s: get_larger(&(r.start), &(i.interval().start)), 
                                        e: get_smaller(&(r.end), &(i.interval().end))});
                n += 1;
            }
            if n >= 1 {
                tmp_v.sort_by(|a, b|
                    if a.s == b.s {
                        a.e.cmp(&(b.e))
                    } else {
                        a.s.cmp(&(b.s))
                    }
                );
                let mut s = tmp_v[0].s;
                let mut e = tmp_v[0].e;
                if n >= 2 {
                    for v in &mut tmp_v[1..(n as usize)] {
                        if v.s < s { s = v.s; }
                        if v.e > e { e = v.e; }
                    }
                }
                let bed = BedRecord::new(chr.clone(), s, e, r.attr.clone());
                new_bed.push(bed);
                tmp_v.clear();
            }
        }
    }
    new_bed
}



#[cfg(test)]
mod tests {
    use crate::intersect::*;

    #[test]
    fn test_intersect() {
        let bed1 = BedRecord::new("chr1".to_string(), 100_u32, 200_u32, vec!["aa", "bb"]);
        let bed2 = BedRecord::new("chr1".to_string(), 150_u32, 250_u32, vec!["cc", "dd"]);
        let bed3 = BedRecord::new("chr2".to_string(), 100_u32, 200_u32, vec!["ccc", "ddd"]);
        let bed4 = BedRecord::new("chr2".to_string(), 150_u32, 250_u32, vec!["ccc", "ddd"]);
        let bed5 = BedRecord::new("chr1".to_string(), 120_u32, 160_u32, [11, 22]);
        let bed6 = BedRecord::new("chr1".to_string(), 110_u32, 145_u32, [33, 44]);
        let bed7 = BedRecord::new("chr2".to_string(), 120_u32, 130_u32, [55, 66]);
        
        let mut bed_a = BedRecords::new();
        let mut bed_b = BedRecords::new();
        bed_a.push(bed1);
        bed_a.push(bed2);
        bed_a.push(bed3);
        bed_a.push(bed4);
        bed_b.push(bed5);
        bed_b.push(bed6);
        bed_b.push(bed7);
        
        // test wa, bed_b should be sorted (otherwise, it will be sorted inside of intersection func)
        bed_b.sort_allchr();  // recommend to sort before intersect
        let (chrs, intersect) = intersect_wa(&mut bed_a, &mut bed_b, INTERSECT_DROP_DUPLICATE);
        for chr in &chrs {
            match intersect.get(chr) {
                Some(v) => {
                    for record in v {
                        println!("{:?}", record);
                    }
                },
                None => println!("{chr} was not found in HashMap")
            }
        }
        println!("");
        
        // test wa wb
        let (chrs, intersect) = intersect_wawb(&mut bed_a, &mut bed_b);
        for chr in &chrs {
            match intersect.get(chr) {
                Some(v) => {
                    for (ra, rb) in v {
                        println!("{:?} : {:?}", ra, rb);
                    }
                },
                None => println!("{chr} was not found in HashMap")
            }
        }
        println!("");
        
        // test v
        let (chrs, intersect) = intersect_v(&mut bed_a, &mut bed_b);
        for chr in &chrs {
            match intersect.get(chr) {
                Some(v) => {
                    for record in v {
                        println!("{:?}", record);
                    }
                },
                None => println!("{chr} was not found in HashMap")
            }
        }
        println!("");
        
        // test merge
        let bed_m = merge(&mut bed_a);
        for chr in &(bed_m.sorted_chrs) {
            match bed_m.map.get(chr) {
                Some(v) => {
                    for record in v {
                        println!("{:?}", record);
                    }
                },
                None => println!("{chr} was not found in HashMap")
            }
        }
        println!("");
    
        // test a
        let mut bed_ra = intersect_a(&mut bed_a, &mut bed_b);
        bed_ra.sort_chrnames();
        for chr in &(bed_ra.sorted_chrs) {
            match bed_ra.map.get(chr) {
                Some(v) => {
                    for record in v {
                        println!("{:?}", record);
                    }
                },
                None => println!("{chr} was not found in HashMap")
            }
        }
        println!("");
    }
} // mod tests

use bio::data_structures::interval_tree::{IntervalTree};
use std::collections::HashMap;
use std::ptr;


/*
 Stores one bed line.
 */
#[derive(Debug)]
pub struct BedRecord<T> {
    chr:    String,
    start:  u64,
    end:    u64,
    attr:   T,
}

impl<T> BedRecord<T> {
    pub fn new(chr: String, start: u64, end: u64, attr: T) -> Self {
        BedRecord{chr: chr, start: start, end: end, attr: attr}
    }
}



/*
 Stores multiple BedRecord objects.
 */
#[derive(Debug)]
pub struct BedRecords<T> {
    map: HashMap<String, Vec<BedRecord<T>>>,  // (name of chr, vec)
    chrs: Vec<String>,
    is_sorted: HashMap<String, bool>,
    sorted_chrs: Vec<String>,
    is_chr_sorted: bool,
}

impl<T> BedRecords<T> {
    pub fn new() -> Self {
        BedRecords{map: HashMap::<String, Vec<BedRecord<T>>>::new(), 
                   chrs: Vec::<String>::new(), 
                   sorted_chrs: Vec::<String>::new(),
                   is_sorted: HashMap::<String, bool>::new(),
                   is_chr_sorted: false}
    }
    
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
    
    pub fn sort_chrnames(&mut self) {
        if ! self.is_chr_sorted {
            self.sorted_chrs = self.chrs.to_vec();
            self.sorted_chrs.sort();
            self.is_chr_sorted = true;
        }
    }
    
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



pub fn intersect_wa<'a, T1, T2>(a: &'a mut BedRecords<T1>, 
                                b: &'a mut BedRecords<T2>, 
                                allow_duplicate: bool) 
        -> (Vec<String>, HashMap::<String, Vec<&'a BedRecord<T1>>>) {
    // make a vector
    let mut intersect_map: HashMap<String, Vec<&'a BedRecord<T1>>> = HashMap::new();
    let mut intersect_chrs: Vec<String> = Vec::new();

    // construct tree
    let mut tree;
    
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
        
        tree = IntervalTree::new();
        for r in &(a.map[chr]) {
            tree.insert(r.start .. r.end, r);
        }
        for r in &(b.map[chr]) {
            for i in tree.find(r.start .. r.end) {
                match intersect_map.get_mut(chr) {
                    Some(v) => {
                        if allow_duplicate {
                            v.push(i.data());
                        } else if ! ptr::eq(*(i.data()), *(v.last().unwrap())) {
                            v.push(i.data());
                        }
                    },
                    None => {
                        let mut v: Vec<&'a BedRecord<T1>> = Vec::new();
                        v.push(i.data());
                        intersect_map.insert(chr.to_string(), v);
                        intersect_chrs.push(chr.to_string());
                    }
                }
            }
        }
    }
    (intersect_chrs, intersect_map)
}


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



#[cfg(test)]
mod tests {
    use crate::intersect::*;

    #[test]
    fn test_intersect() {
        let bed1 = BedRecord::new("chr1".to_string(), 100, 200, vec!["aa", "bb"]);
        let bed2 = BedRecord::new("chr1".to_string(), 150, 250, vec!["cc", "dd"]);
        let bed3 = BedRecord::new("chr2".to_string(), 100, 200, vec!["ccc", "ddd"]);
        let bed4 = BedRecord::new("chr2".to_string(), 150, 250, vec!["ccc", "ddd"]);
        let bed5 = BedRecord::new("chr1".to_string(), 120, 140, vec![11, 22]);
        let bed6 = BedRecord::new("chr2".to_string(), 120, 130, vec![111, 222]);

        let mut bed_a = BedRecords::new();
        let mut bed_b = BedRecords::new();
        bed_a.push(bed1);
        bed_a.push(bed2);
        bed_a.push(bed3);
        bed_a.push(bed4);
        bed_b.push(bed5);
        bed_b.push(bed6);
        
        // test wa, bed_b should be sorted (otherwise, it will be sorted inside of intersection func)
        let ALLOW_DUPLICATE = false;  // default should be false
        bed_b.sort_allchr();  // recommend to sort before intersect
        let (chrs, intersect) = intersect_wa(&mut bed_a, &mut bed_b, ALLOW_DUPLICATE);
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
    }

} // mod tests

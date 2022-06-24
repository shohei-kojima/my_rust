 
#[inline]
fn utf8_to_string(v: &[u8]) -> String {
    String::from_utf8(v.to_vec()).unwrap()
}

#[inline]
fn utf8_to_str(v: &[u8]) -> &str {
    std::str::from_utf8(v).unwrap()
}

pub fn read_bam_header(f: &str) {
    use rust_htslib::bam::Read;
    
    let bam    = rust_htslib::bam::Reader::from_path(f).unwrap();
    
    // let header = rust_htslib::bam::Header::from_template(bam.header());
    // let header_string = String::from_utf8(bam.header().to_bytes()).unwrap();
    // println!("{}", header_string);
    
    // let record = bam.header().as_bytes().to_owned();
    // println!("{}", String::from_utf8(record).unwrap());
    
    let header = rust_htslib::bam::Header::from_template(bam.header());
    for (key, records) in header.to_hashmap() {
        println!("key: {}", key);
        for record in records {
            for (k, v) in record.iter() {
                println!("@{}\t{}:{}", key, k, v);
            }
        }
    }
}

pub fn read_bam(f: &str) {
    use rust_htslib::bam::Read;
        
    let mut bam = rust_htslib::bam::Reader::from_path(f).unwrap();
    bam.set_threads(2).unwrap();
    // bam.set_reference(crai)  // for CRAM
    let mut n = 0;
    for read in bam.rc_records().map(|x| x.expect("Failure parsing Bam file")) {
        // https://docs.rs/rust-htslib/0.19.0/rust_htslib/bam/record/struct.Record.html
        println!("tid: {}: {}", read.tid(), crate::utils::type_of(read.tid()));  // i32
        println!("qname: {}", utf8_to_str(read.qname()));
        println!("start: {}: {}", read.pos(), crate::utils::type_of(read.pos()));  // i32
        println!("flag: {}: {}", read.flags(), crate::utils::type_of(read.flags()));  // u16
        println!();
        n += 1;
        if n == 2 {
            return;
        }
    }
}

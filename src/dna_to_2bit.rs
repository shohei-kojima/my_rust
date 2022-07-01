// Author: Shohei Kojima @ RIKEN

//! Converts DNA to 2bit.


static dna_to_2bitf_u64: [u64; 256] = {
    let shift = 62;
    let mut arr: [u64; 256] = [0; 256];
    arr['A' as u8 as usize] = 0; arr['a' as u8 as usize] = 0;
    arr['T' as u8 as usize] = 1; arr['t' as u8 as usize] = 1;
    arr['G' as u8 as usize] = 2; arr['g' as u8 as usize] = 2;
    arr['C' as u8 as usize] = 3; arr['c' as u8 as usize] = 3;
    arr
};

static dna_to_2bitr_u64: [u64; 256] = {
    let shift = 62;
    let mut arr: [u64; 256] = [0; 256];
    arr['A' as u8 as usize] = 0 << shift; arr['a' as u8 as usize] = 0 << shift;
    arr['T' as u8 as usize] = 1 << shift; arr['t' as u8 as usize] = 1 << shift;
    arr['G' as u8 as usize] = 2 << shift; arr['g' as u8 as usize] = 2 << shift;
    arr['C' as u8 as usize] = 3 << shift; arr['c' as u8 as usize] = 3 << shift;
    arr
};



/// This converts any-nt DNA to 2bit and stores as u64 in the vector. Reads both forward and reverse strands.
pub fn dna_to_2bit_bidirectional_64(seq: &String, kmer_size: &usize) -> Vec<u64> {
    let mut v: Vec<u64> = Vec::new();
    let seqlen = seq.len();
    if seqlen < kmer_size { return v; }
    let bseq = seq.as_bytes();

    // shift bits
    if (2 * kmer_size) > 64 { return v; }
    let shift_bit = 64 - (2 * kmer_size);

    // first window size
    let mut bit2f: u64 = 0;
    let mut bit2r: u64 = 0;
    let mut nn: usize = 0;
    for i in 0 .. kmer_size {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_u64[bseq[i] as usize];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_u64[bseq[i] as usize];
        if bseq[i] == 78 || bseq[i] == 110 {  // ignore when N or n appears; N = 78, n = 110
            nn = kmer_size - 1;
        } else if nn > 0 {  // within window_size-nt from N or n
            nn -= 1;
        }
    }
    bit2r >>= shift_bit;
    if nn == 0 {
        if bit2f == bit2r {
            v.push(bit2f);
        } else {
            v.push(bit2f);
            v.push(bit2r);
        }
    }
    
    // rolling calc.
    for i in kmer_size .. seqlen {
        // shift back
        bit2r <<= shift_bit;
        // usual processing
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_u64[bseq[i] as usize];
        bit2r >>= 2;
        bit2r |= dna_to_2bitr_u64[bseq[i] as usize];
        // shift forward
        bit2f <<= shift_bit;
        bit2f >>= shift_bit;
        bit2r >>= shift_bit;
        if bseq[i] == 78 || bseq[i] == 110 {  // ignore when N or n appears; N = 78, n = 110
            nn = kmer_size - 1;
        } else if nn > 0 {  // within window_size-nt from N or n
            nn -= 1;
        } else if bit2f == bit2r {
            v.push(bit2f);
        } else {
            v.push(bit2f);
            v.push(bit2r);
        }
    }
    
    v
}


/// This converts any-nt DNA to 2bit and stores as u64 in the vector. Only reads forward strand.
pub fn dna_to_2bit_monodirectional_64(seq: &String, kmer_size: &usize) -> Vec<u64> {
    let mut v: Vec<u64> = Vec::new();
    let seqlen = seq.len();
    if seqlen < kmer_size { return v; }
    let bseq = seq.as_bytes();

    // shift bits
    if (2 * kmer_size) > 64 { return v; }
    let shift_bit = 64 - (2 * kmer_size);

    // first window size
    let mut bit2f: u64 = 0;
    let mut nn: usize = 0;
    for i in 0 .. kmer_size {
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_u64[bseq[i] as usize];
        if bseq[i] == 78 || bseq[i] == 110 {  // ignore when N or n appears; N = 78, n = 110
            nn = kmer_size - 1;
        } else if nn > 0 {  // within window_size-nt from N or n
            nn -= 1;
        }
    }
    if nn == 0 {
        v.push(bit2f);
    }
    
    // rolling calc.
    for i in kmer_size .. seqlen {
        // usual processing
        bit2f <<= 2;
        bit2f |= dna_to_2bitf_u64[bseq[i] as usize];
        // shift forward
        bit2f <<= shift_bit;
        bit2f >>= shift_bit;
        if bseq[i] == 78 || bseq[i] == 110 {  // ignore when N or n appears; N = 78, n = 110
            nn = kmer_size - 1;
        } else if nn > 0 {  // within window_size-nt from N or n
            nn -= 1;
        } else {
            v.push(bit2f);
        }
    }
    
    v
}



/// Implement iterator that takes a seq and returns 2bit k-mers.
// pub struct BitSeqF<'a> {
//     curr: u64,
//     next: u64,
//     bytes: &'a [u8],
//     pos: usize,
//     nn: usize,
// }
// 
// impl<'a> BitSeqF<'a> {
//     pub fn new(seq: &'a str, kmer_size: &usize) -> Option<Self> {
//         let seqlen = seq.len();
//         if seqlen < kmer_size { return None; }
//         let bseq = seq.as_bytes();
// 
//         // shift bits
//         if (2 * kmer_size) > 64 { return None; }
//         let shift_bit = 64 - (2 * kmer_size);
// 
//         // first window size
//         let mut bit2f: u64 = 0;
//         let mut nn: usize = 0;
//         for i in 0 .. kmer_size {
//             bit2f <<= 2;
//             bit2f |= dna_to_2bitf_u64[bseq[i] as usize];
//             if bseq[i] == 78 || bseq[i] == 110 {  // ignore when N or n appears; N = 78, n = 110
//                 nn = kmer_size - 1;
//             } else if nn > 0 {  // within window_size-nt from N or n
//                 nn -= 1;
//             }
//         }
//         Some( BitSeqF { curr: bit2f, next: 0, bytes: bseq, pos: kmer_size, nn: nn } )
//     }
// }





#[cfg(test)]
mod tests {
    use crate::dna_to_2bit::*;
    
    // #[test]
    fn test_2bit() {
        let dna = "ACGATCGACTAACGGATCGAGGCGGCGATCATTTCGAGCTAGGACAACATCACACGCGCGATCGAT".to_string();
        let v = dna_to_2bit_bidirectional_64(&dna, 32);
        println!("{:?}", v.len());
        
        let v = dna_to_2bit_monodirectional_64(&dna, 32);
        println!("{:?}", v.len());
    }
} // mod tests

// Author: Shohei Kojima @ RIKEN

//! DNA sequence manipulations.

/// To convert from lower case DNA/RNA to upper case
pub static TO_UPPER_ARR: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 65,
    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
    90, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
    141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
    179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
    198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,
    217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235,
    236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
    255,
];

/// To convert from upper case DNA/RNA to lower case
pub static TO_LOWER_ARR: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 97, 98, 99, 100, 101, 102, 103,
    104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
    91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
    131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
    169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187,
    188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206,
    207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
    226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244,
    245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
];

/// To convert DNA to complementary (will not do reverse!) e.g. AT-N-gc -> TA-N-cg
pub static TO_COMPLEMENT_ARR: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
    116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
    97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
    173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
    211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
    230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248,
    249, 250, 251, 252, 253, 254, 255,
];

pub fn to_upper_seq(seq: &mut String) {
    unsafe {
        for i in 0..seq.len() {
            seq.as_mut_vec()[i] = TO_UPPER_ARR[seq.as_mut_vec()[i] as usize];
        }
    }
}

pub fn to_lower_seq(seq: &mut String) {
    unsafe {
        for i in 0..seq.len() {
            seq.as_mut_vec()[i] = TO_LOWER_ARR[seq.as_mut_vec()[i] as usize];
        }
    }
}

pub fn to_comp_seq(seq: &mut String) {
    unsafe {
        for i in 0..seq.len() {
            seq.as_mut_vec()[i] = TO_COMPLEMENT_ARR[seq.as_mut_vec()[i] as usize];
        }
    }
}

pub fn to_rev_seq(seq: &mut String) {
    unsafe {
        seq.as_mut_vec().reverse();
    }
}

pub fn to_rev_comp_seq(seq: &mut String) {
    to_rev_seq(seq);
    to_comp_seq(seq);
}

#[cfg(test)]
mod tests {
    use crate::nucl::*;

    // #[test]
    fn test_nucl() {
        let dna = "ATCGAC--Nnn--actagcat".to_string();

        let mut tmp = dna.clone();
        to_upper_seq(&mut tmp);
        println!("{} -> upper: {}", dna, tmp);

        let mut tmp = dna.clone();
        to_lower_seq(&mut tmp);
        println!("{} -> lower: {}", dna, tmp);

        let mut tmp = dna.clone();
        to_comp_seq(&mut tmp);
        println!("{} -> comp: {}", dna, tmp);

        let mut tmp = dna.clone();
        to_rev_seq(&mut tmp);
        println!("{} -> rev: {}", dna, tmp);

        let mut tmp = dna.clone();
        to_rev_comp_seq(&mut tmp);
        println!("{} -> rev_comp: {}", dna, tmp);
    }
} // mod tests

// Just a memo
// fn make_to_complement_arr() -> [u8; 256] {
//     let mut to_comp_arr: [u8; 256] = [0; 256];
//     for n in 0 .. 256 {
//         to_comp_arr[n as usize] = n as u8;
//     }
//     to_comp_arr[65] = 84;  // A -> T
//     to_comp_arr[84] = 65;  // T -> A
//     to_comp_arr[71] = 67;  // G -> C
//     to_comp_arr[67] = 71;  // C -> G
//     to_comp_arr[97] = 116;  // a -> t
//     to_comp_arr[116] = 97;  // t -> a
//     to_comp_arr[103] = 99;  // g -> c
//     to_comp_arr[99] = 103;  // c -> g
//     to_comp_arr
// }

// Just a memo
// fn make_to_upper_arr() -> [u8; 256] {
//     let mut to_upper_arr: [u8; 256] = [0; 256];
//     for n in 0 .. 256 {
//         to_upper_arr[n as usize] = n as u8;
//     }
//     for lower in 97 .. 123 {
//         to_upper_arr[lower as usize] = (lower - 32) as u8;
//     }
//     to_upper_arr
// }

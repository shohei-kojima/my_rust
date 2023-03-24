// Author: Shohei Kojima @ RIKEN

//! Misc functions.

use std::fs;
use std::path::PathBuf;

pub fn type_of<T>(_: T) -> String {
    std::any::type_name::<T>().to_string()
}

pub fn path_exists(fname: &str) -> bool {
    PathBuf::from(fname).exists()
}

pub fn is_file(fname: &str) -> bool {
    PathBuf::from(fname).is_file()
}

pub fn is_dir(fname: &str) -> bool {
    PathBuf::from(fname).is_dir()
}

pub fn makedirs_exist_ok(path: &str) -> std::io::Result<()> {
    if !path_exists(path) {
        fs::create_dir_all(PathBuf::from(path))?;
    } else if !is_dir(path) {
        panic!("{path} is not directory.");
    }
    Ok(())
}

/// Check whether a file exists.
/// Exit with 1 if the file does not exist.
pub fn file_exist_check(file_path: &str) {
    if !path_exists(file_path) {
        panic!("Input file ({}) does not exist.", file_path);
    } else {
        println!("Input file ({}) found.", file_path);
    }
}

/// Check whether a file does not exist.
pub fn file_absence_check(file_path: &str, overwrite: u8) {
    if path_exists(file_path) {
        if overwrite == 0 {
            println!("Warn: Output will be overwritten in {}.", file_path);
        } else {
            panic!("Output file ({}) already exists.", file_path);
        }
    } else {
        println!("Output will be written in ({}).", file_path);
    }
}

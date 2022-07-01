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
    if ! path_exists(path) {
        fs::create_dir_all(PathBuf::from(path))?;
    } else if ! is_dir(path) {
        panic!("{path} is not directory.");
    }
    Ok(())
}

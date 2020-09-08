//! This is a crate for python read until api.
#![crate_name = "dyss"]
#![crate_type = "lib"]
extern crate bio;
extern crate dtw;
extern crate libc;
extern crate rand;
extern crate squiggler;
extern crate utility;
extern crate rayon;
pub mod dyss;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2+2,4);
    }
}


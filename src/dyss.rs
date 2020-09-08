use libc;
use std;
use dtw;
use rayon::prelude::*;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use squiggler;
use utility::utilities;
const WINDOWS:[usize;2] = [4, 8];
const THRESHOLDS:[f64;2] =[1.5, 9.0];
const DELTA:f64 = 0.2;
#[derive(Debug)]
#[repr(C)]
pub struct Dyss {
    threshold: f32,
    num_scouts: usize,
    num_packs: usize,
    querysize: usize,
    power: f32,
    reference: Vec<f32>,
}

impl Dyss {
    fn new(
        num_scouts: usize,
        num_packs: usize,
        ref_path: &str,
        model_path: &str,
        param_path: &str,
        power: u32,
        querysize: usize,
        refsize: usize,
    ) -> Option<Self> {
        let (threshold, _specificity): (f32, f32) = match Self::get_threshold(
            &Path::new(param_path),
            refsize,
            power,
            num_packs,
            num_scouts,
        ) {
            Some(res) => res,
            None => return None,
        };
        eprintln!("{},{},{},{},{}", refsize, power, num_packs, num_scouts,threshold);
        let model = match squiggler::Squiggler::new(&Path::new(model_path)) {
            Ok(res) => res,
            Err(_why) => {
                eprintln!(r#"
Something went wrong during the conversion from fasta file into raw-signal-like pattern.
First of all, please check the required files(model file and fasta file) exists.
Next, please make sure that both ont_fast5_api and h5py library have been installed correctly into your PYTHON3.
It might be because the model file is not written in a valid format. Please use ont_kmer_model as is.
Or feel free to contact us: banmasutani@gmail.com
"#);
                return None;
            }
        };
        let power = power as f32 / 100.;
        let (temp,rev) = match utilities::setup_template_complement_autofit(&Path::new(ref_path),refsize){
            Ok((temp,rev)) => (temp,rev),
            Err(why) => {
                eprintln!("{:?}",why);
                return None;
            }
        };
        let temp = utilities::convert_to_squiggle(&temp,&model);
        let rev = utilities::convert_to_squiggle(&rev,&model);
        let reference:Vec<_> = temp.into_iter().chain(rev.into_iter()).collect();
        eprintln!("reflen:{:?}",reference.len());
        Some(Dyss {
            threshold,
            num_scouts,
            num_packs,
            reference,
            querysize,
            power,
        })
    }
    #[inline]
    fn get_threshold(
        param_path: &Path,
        refsize: usize,
        power: u32,
        num_packs: usize,
        num_scouts: usize,
    ) -> Option<(f32, f32)> {
        match File::open(param_path).and_then(|e| {
            BufReader::new(e).lines()
                .filter_map(|e|e.ok())
                .skip(1)//header
                .filter_map(|e| Self::parse_line(e))
                .filter(|e| e.0 == refsize && e.1 == power &&
                        e.2 == num_packs && e.3 == num_scouts)
                .nth(0)
                .ok_or(std::io::Error::from(std::io::ErrorKind::Other))
        }) {
            Ok(res) => Some((res.4, res.5)),
            Err(why) => {
                eprintln!("{:?}", why);
                None
            }
        }
    }
    #[inline]
    fn parse_line(line: String) -> Option<(usize, u32, usize, usize, f32, f32)> {
        let contents: Vec<_> = line.split(',').collect();
        let refsize: usize = match contents[0].parse::<usize>() {
            Ok(res) => res * 1000,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        let power: u32 = match contents[1].parse() {
            Ok(res) => res,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        let num_packs: usize = match contents[2].parse() {
            Ok(res) => res,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        let num_scouts: usize = match contents[3].parse() {
            Ok(res) => res,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        let threshold: f32 = match contents[4].parse() {
            Ok(res) => res,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        let specificity: f32 = match contents[5].parse() {
            Ok(res) => res,
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        };
        Some((
            refsize,
            power,
            num_packs,
            num_scouts,
            threshold,
            specificity,
        ))
    }
    #[inline]
    fn classify(&self, query: &[i32]) -> i32 {
        use squiggler;
        let before = query.len();
        let query = utilities::trim_front_vec(&query,70,15);
        let query = utilities::event_detect_mult_ttest(&query,&WINDOWS,&THRESHOLDS,DELTA);
        if query.len() < self.querysize {
            return 2;// not enough length
        }

        let query:Vec<_> = query.into_iter().map(|e|e.mean as f32).collect();
        let query = dtw::normalize(&query, dtw::NormalizeType::Z);
        let query = squiggler::dedup(&query, self.power);
        if query.len() < self.querysize {
            eprintln!("query chunked {},{}",before,query.len());
            return 2;// not enough length
        }
        // let preprocess = time::Instant::now();
        if let Ok((score, _, _)) = dtw::scouting_threshold_dtw(
            &query[0..self.querysize],
            &self.reference,
            &hill,
            Some(self.num_scouts),
            Some(self.num_packs),
            self.threshold
        ){
            // let process = time::Instant::now();
            // eprintln!("OK {:?},{:?}",preprocess-start ,process - preprocess);
            eprintln!("score,{}", score);
            if score < self.threshold {
                1
            }else{
                0
            }
        } else {
            // let process = time::Instant::now();
            // eprintln!("ER {:?},{:?}",preprocess-start ,process - preprocess);
            0
        }
    }
}

#[inline]
fn hill(x: &f32, y: &f32) -> f32 {
    let d = (x - y).powi(2);
    d / (1.0 + d)
}

#[no_mangle]
pub extern "C" fn construct_dyss(
    num_scouts: *const libc::c_int,
    num_packs: *const libc::c_int,
    ref_path: *const libc::c_char,
    model_path: *const libc::c_char,
    param_path: *const libc::c_char,
    power: libc::c_int,
    querysize: libc::size_t,
    refsize: libc::size_t,
) -> *mut Dyss {
    use std::ffi::CStr;
    let ref_path = unsafe {
        match CStr::from_ptr(ref_path).to_str() {
            Ok(res) => res,
            Err(_why) => {
                eprintln!(r#"
An invalid reference file was sapplied.
DySS halted its execution. Please check your command line arguments(For example, confirm the file exists).
Or contact us: banmasutani@gmail.com
"#);
                return std::ptr::null_mut();
            }
        }
    };
    let model_path = unsafe {
        match CStr::from_ptr(model_path).to_str() {
            Ok(res) => res,
            Err(_why) => {
                eprintln!(r#"
An invalid model file was sapplied.
DySS halted its execution. Please check your command line arguments(For example, confirm the file exists).
Or contact us: banmasutani@gmail.com
"#);
                return std::ptr::null_mut();
            }
        }
    };
    let param_path = unsafe {
        match CStr::from_ptr(param_path).to_str() {
            Ok(res) => res,
            Err(_why) => {
                eprintln!(r#"
An invalid model file was sapplied.
DySS halted its execution. Please check your command line arguments(For example, confirm the file exists).
Or contact us: banmasutani@gmail.com
"#);
                return std::ptr::null_mut();
            }
        }
    };
    let res = match Dyss::new(
        num_scouts as usize,
        num_packs as usize,
        ref_path,
        model_path,
        param_path,
        power as u32,
        querysize as usize,
        refsize as usize,
    ) {
        Some(res) => res,
        None => {
            eprintln!(r#"
DySS could not be constructed properly. Execution halt.
Please check whether the files exist.
If the problem persists, plese contact me: banmasutani@gmail.com
"#);
            return std::ptr::null_mut();
        }
    };
    let res = Box::into_raw(Box::new(res));
    res
}

#[no_mangle]
pub extern "C" fn dyss_classify(
    ptr: *const Dyss,
    query: *const libc::c_int,
    length: libc::size_t,
) -> i32 {
    //    let start = Instant::now();
    let query = if query.is_null() {
        //eprintln!("query is null");
        return 0;
    } else {
        unsafe{std::slice::from_raw_parts(query, length as usize) }
    };
    //    eprint!("querysetup:{}",start.elapsed().subsec_nanos()/(10u32.pow(3)));
    let classifier = if ptr.is_null() {
        //eprintln!("classifier is null. Maybe you free this struct previously");
        return 0;
    } else {
        unsafe { & *ptr }
    };
    //    let c_start = Instant::now();
    classifier.classify(& query)
}

#[no_mangle]
pub extern "C" fn dyss_destructor(ptr: *mut Dyss) {
    let _classifier = if ptr.is_null() {
        eprintln!("classifier is null. Maybe you free this struct previously");
        return;
    } else {
        unsafe { Box::from_raw(ptr) }
    };
}

#[no_mangle]
pub extern "C" fn is_null(ptr: *mut Dyss) -> i32 {
    if ptr.is_null() {
        0
    } else {
        1
    }
}

#[test]
fn parameter_parse() {
    let (threshold, specificity) = Dyss::get_threshold(
        &Path::new("/home/ban-m/work/irabu/data/parameters.csv"),
        100_000,
        20,
        1,
        10,
    ).unwrap();
    assert!(threshold > 21.6);
}

#[no_mangle]
pub extern "C" fn batch_classify(dyss:*const Dyss,data:*const *const libc::c_int,lengths:*const libc::size_t,
                                 num_of_data:libc::size_t,result:*mut libc::c_int)->libc::c_int{
    let dyss = if dyss.is_null(){
        return 0
    }else{
        unsafe{& *dyss}
    };
    let lengths = if lengths.is_null(){
        return 2
    }else{
        unsafe{std::slice::from_raw_parts(lengths,num_of_data)}
    };
    let data = if data.is_null(){
        return 3
    }else{
        unsafe{std::slice::from_raw_parts(data,num_of_data)}
    };
    let mut chunks = Vec::with_capacity(num_of_data);
    for index in 0..num_of_data{
        if data[index].is_null(){
            return 4
        }else{
            chunks.push(unsafe{std::slice::from_raw_parts(data[index],lengths[index])});
        }
    }
    let result = if result.is_null(){
        return 5
    }else{
        unsafe{std::slice::from_raw_parts_mut(result,num_of_data)}
    };
    let classified:Vec<_> = chunks.par_iter().map(|&e|dyss.classify(e)).collect();
    for index in 0..num_of_data{
        result[index] = classified[index];
    }
    1
}

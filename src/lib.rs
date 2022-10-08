#![allow(unused)]

/// We use Levenshtein distance on 2-combination
/// of all the sequences after running them through
/// n_trim().
/// Then we filter the distances that are LOWER than
/// the threshold and we KEEP them.
/// After that we use fasthash to get their u64 hash
/// and if they are equal, we keep them.
/// We flatten the iterator first.
/// Then we collect the string vector.
/// Contact Chubak#7400 if you got any questions.

#[macro_use]
extern crate lazy_static;

use fasthash::metro;
use lazy_static::__Deref;
use pyo3::prelude::*;
use hashbrown::HashSet;
use rayon::prelude::*;
use std::sync::Arc;
use parking_lot::Mutex;
use std::ops::DerefMut;
use bio::alphabets::rna::revcomp;
use std::fs::read_to_string;
use std::time::{SystemTime, Duration};
use std::thread;
use crossbeam_channel::{Sender, unbounded};
use waker_fn::waker_fn;
use threadpool::ThreadPool;

// these functions are string distancde and the old revcom

/* 
#[derive(Clone)]
struct SiftTrans {
    pub c_a: isize,
    pub c_b: isize,
    pub is_trans: bool,
}

fn sift_four_common(a: &String, b: &String, max_offset: isize, max_distance: Option<f32>) -> f32 {
 
    if a.len() == 0 {
        if b.len() == 0 {
            return 0f32;
        }

        return b.len() as f32;
    }

    let len_a = a.len() as isize;
    let len_b = b.len() as isize;

    let mut c_a = 0isize;
    let mut c_b = 0isize;

    let mut lcss = 0;
    let mut local_cs = 0;
    let mut trans = 0;

    let mut debug = 0;

    let mut offset_vec: Vec<SiftTrans> = vec![];

    while (c_a < len_a) && (c_b < len_b) {     
        let char_a = a.chars().nth(c_a as usize).unwrap() as u8;
        let char_b = b.chars().nth(c_b as usize).unwrap() as u8;

        if char_a == char_b {
            local_cs += 1;

            let mut is_trans = false;

            offset_vec
                .clone()
                .iter_mut()
                .enumerate()
                .for_each(|(i, ofs)| {
                    if c_a <= ofs.c_a || c_b <= ofs.c_b {
                        is_trans = (c_b - c_a).abs() >= (ofs.c_b - ofs.c_a).abs();

                        if is_trans {
                            trans += 1;
                        } else {
                            if !ofs.is_trans {
                                ofs.is_trans = true;
                                trans += 1;
                            }
                        }
                    } else {
                        if (c_a > ofs.c_b && c_b > ofs.c_a) {
                            offset_vec.splice(i..i + 1, vec![]);
                        } else {
                            trans += 1;
                        }
                    }
                });

            offset_vec.push(SiftTrans { c_a, c_b, is_trans });
        } else {
            lcss += local_cs;
            local_cs = 0;

            if c_a != c_b {
                c_b = isize::min(c_a, c_b);
                c_a = c_b;
            }

            let mut i = 0isize;

            loop {
                if i < max_offset && (c_a + i < len_a || c_b + i < len_b) {
                    break;
                }

                if (c_a + i < len_a)
                    && (a.chars().nth((c_a + i) as usize) == b.chars().nth(c_b as usize))
                {
                    c_a += i - 1;
                    c_b -= 1;

                    break;
                }

                if (c_b + i < len_b)
                    && (a.chars().nth(c_a as usize) == b.chars().nth((c_b + i) as usize))
                {
                    c_a -= 1;
                    c_b += i - 1;
                    break;
                }

                i += 1;
            }
        }

        c_a += 1;
        c_b += 1;

        if max_distance.is_some() {
            let max_dist = max_distance.unwrap();

            let temporary_distance = (isize::max(c_a, c_b) as f32) - (lcss as f32) + (trans as f32);

            if temporary_distance >= max_dist {
                return temporary_distance.round();
            }
        }

        if (c_a >= len_a) || (c_b >= len_b) {
            lcss += local_cs;
            local_cs = 0;

            c_b = isize::min(c_a, c_b);
            c_a = c_b;
        }
    }

    lcss += local_cs;

    let ret = (isize::max(len_a, len_b) - lcss + trans) as f32;

    ret.round()
}


fn reverse_complement(line: &String) -> String {
    String::from_iter(
        line.chars()
            .map(|ch| {
                let ret = {
                    if vec!['A', 'T', 'G', 'C'].contains(&ch) {
                        return *COMPS.get(&ch).unwrap();
                    }
                    ch
                };
                ret
            })
            .collect::<Vec<char>>(),
    )
}


fn reverse_complement(line: &mut String) {
    *line = line.replace("A", "T");
    *line = line.replace("T", "A");
    *line = line.replace("C", "G");
    *line = line.replace("G", "C");    
}


fn n_trim(parent_seq: String, min_seq_length: usize) -> Vec<String> {
    if !parent_seq.chars().any(|x| x == 'N') {
        return vec![parent_seq];
    }

    let mut filtered = parent_seq
        .char_indices()
        .filter(|(i, ch)| *ch == 'N')
        .map(|(i, _)| i)
        .collect::<Vec<usize>>();

    filtered.splice(..0, vec![0usize]);
    filtered.push(parent_seq.len() - 1);

    let ret = (0usize..filtered.len() - 1)
        .step_by(2)
        .map(|(i)| {
            let (start, end) = (filtered.get(i).unwrap(), filtered.get(i + 1).unwrap());

            let ret = if (end - start) + 1 >= min_seq_length {
                let mut str = parent_seq[*start..*end].to_string();
                str = str.trim_matches('N').to_string();

                str
            } else {
                return String::new();
            };

            ret
        })
        .filter(|x| *x != "")
        .collect::<Vec<String>>();

    match ret.len() == 0 {
        true => vec![parent_seq],
        false => ret,
    }
}



fn add_header(seq: String, i: &mut usize) -> String {
    *i += 1;

    format!(">NODE_{}_length_{}\n{}\n", i, seq.len(), seq)
}

fn add_revcomp_to_hashset(v: Vec<u8>, sender: Sender<Vec<u8>>) {                        
    
        
}

#[pyfunction]
fn dedup_lines(lines: Vec<String>, min_seq_length: usize) -> Vec<String> {
    let mut hs = DashSet::<Vec<u8>>::new();
    let mut ret = vec![];
    let (s, r) = unbounded();

    (0..lines.len())
        .step_by(2)
        .collect::<Vec<usize>>()
        .into_iter()
        .map(|i| lines.get(i + 1).unwrap().as_bytes().to_vec())
        .for_each(|v|  {
            let ss = s.clone();
            let a_comp = v.clone();

            thread::spawn(move || {
                ss.send(a_comp.clone());
                let rvcmp = revcomp(a_comp);
                ss.send(rvcmp);
            });            
        });

    
    match r.recv() {
        Ok(vv) => {
            hs.insert(vv);
        },
        Err(e) => panic!("{}", e),
    }
    
    
    let (ss, rs) = unbounded();
     
    hs
        .into_iter()
        .enumerate()
        .for_each(|(i, v)| {
            let ssc = ss.clone();

            thread::spawn(move || {
                let st = String::from_utf8(v).unwrap();

                ssc.send(st);
            });            
        });

    let mut index = 0usize;
    
    match rs.recv() {
        Ok(s) => add_header(
            s,
            &mut index,
        ),
        Err(e) => panic!("{}", e),
    };
            
    

    println!("Main dedup operation done, returning strings...");

    ret
}


*/

fn add_header(seq: String, i: usize) -> String {
    format!(">NODE_{}_length_{}\n{}\n", i, seq.len(), seq)
}

#[pyfunction]
fn dedup_lines(lines: Vec<String>, min_seq_length: usize) -> Vec<String> {
    let mut hashset = HashSet::<String>::new();
    let len_vec = (0..lines.len()).step_by(2).into_iter().collect::<Vec<usize>>();

    let n_workers = 16;
    let pool = ThreadPool::new(n_workers);

    let am_hs = Arc::new(Mutex::new(hashset));

    for i in len_vec {
        let s = lines.get(i + 1).unwrap().clone();
        let s_am = Arc::new(Mutex::new(s)); 
        let hs_am = Arc::clone(&am_hs);

        pool.execute(move || {
            let s_mutex_gaurd = s_am.lock();
            let ss = s_mutex_gaurd.deref();

            let mut hs_lock = hs_am.lock();
            let mut hs_deref = hs_lock.deref_mut();
            
            hs_deref.insert(ss.clone());
            
            let rcomp = revcomp(ss.clone().as_bytes().to_vec());
            let rcomp_str = String::from_utf8(rcomp).unwrap();
            hs_deref.insert(rcomp_str);
        })
    }

    println!("Main job done.");

    let mut hs_mg = am_hs.lock();
    let hs_deref = hs_mg.deref();

    let ret = hs_deref
        .into_par_iter()   
        .map(|s| s.clone())
        .collect::<Vec<String>>();

    ret
}



#[pymodule]
fn phymmr_dedup(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dedup_lines, m)?)?;
    Ok(())
}



#[test]
fn test_stif_four() {
    let str = read_to_string("SRR6824218.fa").unwrap();
    let lines = str.lines().into_iter().map(|x| x.to_string()).collect::<Vec<String>>();

    println!("Starting main op");
    let t = SystemTime::now();
    dedup_lines(lines, 90);
    assert_eq!(Duration::from_secs(1), t.elapsed().unwrap());
}   

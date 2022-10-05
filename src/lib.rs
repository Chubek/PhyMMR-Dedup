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
use itertools::Itertools;
use levenshtein::levenshtein;
use pyo3::prelude::*;
use std::collections::HashMap;

lazy_static! {
    static ref COMPS: HashMap<char, char> = {
        let mut comps = HashMap::<char, char>::new();

        comps.insert('A', 'T');
        comps.insert('C', 'G');
        comps.insert('G', 'C');
        comps.insert('T', 'A');

        comps
    };
}

fn reverse_complement(line: String) -> String {
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

fn n_trim(parent_seq: String, min_seq_length: usize) -> Vec<String> {
    if !parent_seq.chars().any(|x| x == 'N') {
        return vec![parent_seq];
    }

    parent_seq
        .char_indices()
        .filter(|(i, ch)| *ch == 'N' || *i == 0 || *i == parent_seq.len() - 2)
        .map(|(i, _)| {
            let (start, end) = (i, i + 1);

            let ret = if (end - start) + 1 >= min_seq_length {
                let mut str = parent_seq[start..end].to_string();
                str = str.trim_matches('N').to_string();

                str
            } else {
                return String::new();
            };

            ret
        })
        .collect::<Vec<String>>()
}

#[pyfunction]
fn dedup_lines(lines: Vec<String>, min_seq_length: usize, leven_thresh: usize) -> Vec<String> {
    let mut header_counter = 1;


    lines
        .into_iter()
        .map(|x| {
            let trim = x.trim().to_string();

            trim
        })
        .filter(|x| !(x.chars().next().unwrap() == '>') && x.len() >= min_seq_length)
        .map(|x| {
            let combine_leven = {
                n_trim(x, min_seq_length)
                    .into_iter()
                    .combinations(2)
                    .into_iter()
                    .map(|v| {
                        let (a, b) = (v.get(0).unwrap(), v.get(1).unwrap());
                        
                        let a_comp = reverse_complement(a.clone());
                        let b_comp = reverse_complement(b.clone());
                        
                        let leven = levenshtein(a, b);
                        let leven_comp_a = levenshtein(&a_comp, b);
                        let leven_comp_b = levenshtein(a, &b_comp);


                        (leven, leven_comp_a, leven_comp_b, v)
                    })
                    .filter(|(
                        leven,
                        leven_comp_a, 
                        leven_comp_b,  
                        _)| *leven <= leven_thresh 
                            && *leven_comp_a <= leven_thresh 
                            && *leven_comp_b <= leven_thresh
                        )
                    .map(|(_, _, _, v)| v)
                    .collect::<Vec<Vec<String>>>()
            };

            combine_leven
        })
        .map(|v| {
            let hashes = {
                v.into_iter()
                    .map(|vi| {
                        let (a, b) = (vi.get(0).unwrap(), vi.get(1).unwrap());

                        let hash_a = metro::hash64(a);
                        let hash_b = metro::hash64(b);

                        let length = a.chars().count();

                        (hash_a, hash_b, a.clone(), length)
                    })
                    .collect::<Vec<(u64, u64, String, usize)>>()
            };

            hashes
        })
        .flatten()
        .filter(|(h1, h2, _, _)| *h1 != *h2)
        .map(|(_, _, seq, length)| {
            let header = format!(">NODE_{}_length_{}", header_counter, length);
            header_counter += 1;
            format!("{}\n'{}'\n", header, seq)
        })
        .collect::<Vec<String>>()
}

#[pymodule]
fn phymmr_dedup(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dedup_lines, m)?)?;
    Ok(())
}

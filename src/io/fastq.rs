use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use io::seq::Seq;
use io::seq::Cleanable;

use regex::Regex;

#[test]
/// Test whether we can read the test data carefully.
fn test_the_careful_reader () {
    use std::fs::File;
    let my_file = File::open("testdata/four_reads.fastq").expect("Could not open file");
    let my_buffer=BufReader::new(my_file);
    let mut fastq_reader=FastqReader::new_careful(my_buffer);
    let seq_obj = fastq_reader.next().expect("Could not read four_reads.fastq");
    assert_eq!(seq_obj.seq.trim(), "AAAGTGCTCTTAACTTGTCCCGCTCCACATCAGCGCGACATCAATCGACATTAAACCGAGTATCTTGTCAGCCTGGGGTGACGATGCGTCCCATTAAAGT");
}

#[test]
/// Test whether we can read the test data quickly.
fn test_the_fast_reader () {
    use std::fs::File;
    let my_file = File::open("testdata/four_reads.fastq").expect("Could not open file");
    let my_buffer=BufReader::new(my_file);
    let mut fastq_reader=FastqReader::new(my_buffer);
    let seq_obj = fastq_reader.next().expect("Could not read four_reads.fastq");
    assert_eq!(seq_obj.seq.trim(), "AAAGTGCTCTTAACTTGTCCCGCTCCACATCAGCGCGACATCAATCGACATTAAACCGAGTATCTTGTCAGCCTGGGGTGACGATGCGTCCCATTAAAGT");
}

/// A FastQ reader
pub struct FastqReader<R: io::Read> {
  reader:           io::BufReader<R>,
  quickly:          bool,
}

impl<R: io::Read> FastqReader<R>{
  pub fn new(reader: R) -> FastqReader<R> {
    FastqReader {
        reader : BufReader::new(reader),
        quickly: true,
    }
  }
  pub fn new_careful(reader: R) -> FastqReader<R> {
    FastqReader {
        reader : BufReader::new(reader),
        quickly: false,
    }
  }
}

impl<R: Read> Iterator for FastqReader<R> {
    type Item = Seq;

    /// There are two flavors of next: either read 
    /// quickly or read carefully, depending on the 
    /// 'quickly' bool.
    fn next(&mut self) -> Option<Seq> {
        // Read a fastq entry but assume that there are only
        // four lines per entry (id, seq, plus, qual).
        if self.quickly {
            let mut id=    String::new();
            let mut seq=   String::new();
            let mut qual=  String::new();

            // Read the ID of the entry
            match self.reader.read_line(&mut id) {
                // if we're expecting an ID line, but
                // there are zero bytes read, then we are
                // at the end of the file. Break.
                Ok(0) => return None,
                Ok(_n) => {},
                Err(error) => {
                    panic!("ERROR: could not read the ID line: {}",error);
                }
            };

            self.reader.read_line(&mut seq).expect("ERROR: could not read sequence line");

            // skip the plus sign
            self.reader.skip_until(b'\n').expect("ERROR: plus sign line not found");

            self.reader.read_line(&mut qual).expect("ERROR: could not read qual line");

            let seq = Seq::new(&id, &seq, &qual);

            return Some(seq)
        }
        // Read a fastq entry in the most correct way possible,
        // allowing for whitespace in seq and qual lines.
        else {
            let whitespace_regex = Regex::new(r"(\s+)").expect("malformed regex");

            let mut id=    String::new();
            let mut seq=   String::new();
            let mut qual=  String::new();

            // Read the ID of the entry
            match self.reader.read_line(&mut id) {
                Ok(0) => return None,
                Ok(_n) => {},
                Err(error) => {
                    panic!("ERROR: {}",error);
                }
              
            }
            // Read the DNA line of the entry and count
            // how long it is.
            'dna: loop{
                // clears the buffer without changing capacity
                match self.reader.read_line(&mut seq) {
                    Ok(0) => panic!("ERROR: incomplete entry (no seq line), seqid {}\nbuf {}", id.trim(), seq),
                    Ok(n) => {
                        // if we hit the qual line, then it is a single
                        // character, +
                        let check_for_plus = seq.len() - n;
                        
                        if &seq[check_for_plus..check_for_plus+1] == "+" {
                            seq.truncate(check_for_plus);
                            break 'dna;
                        }
                    }
                    Err(error) => {
                        panic!("ERROR while reading seq for ID {}: {}",id.trim(),error);
                    }
                }
            }
            // remove all whitespace
            seq = whitespace_regex.replace_all(&seq,"").into_owned();
            let read_length :usize=seq.len(); 
            
            // build onto the qual line until it has the right
            // number of bytes.
            'qual: loop{
                match self.reader.read_line(&mut qual) {
                    Ok(0) => panic!("ERROR: incomplete entry (no qual line), seqid {}\nqual {}", id.trim(),qual),
                    Ok(_n) => {
                        qual = whitespace_regex.replace_all(&qual,"").into_owned();
                        if qual.len() >= read_length {
                          break 'qual;
                        }
                    }
                    Err(error) => {
                        panic!("ERROR while reading qual for ID {}: {}",qual.trim(),error);
                    }
                }
            }

            let seq = Seq::new(&id, &seq, &qual);
            return Some(seq)
        }
    }

}


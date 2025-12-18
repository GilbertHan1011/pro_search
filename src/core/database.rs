use std::fs::File;
use std::io::{BufRead, BufReader,Read};
use std::path::Path;
use anyhow::{Context, Result};

use flate2::read::GzDecoder; // gzip decoder
use zstd::stream::read::Decoder as ZstdDecoder; // zstd decoder

pub struct Database {
    pub ids: Vec<String>,
    pub data: Vec<u8>,        
    pub offsets: Vec<usize>,  // Flattened offsets for quick access
}

impl Database {
    pub fn load_from_fasta<P: AsRef<Path>>(path:P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open file: {:?}", path))?;
        let reader :Box<dyn Read> = match path.extension().and_then(|s| s.to_str()){
            Some("gz") => {
                println!("Detected gzip file");
                Box::new(GzDecoder::new(file))
            }
            Some("zst") => {
                println!("Detected zstd file");
                Box::new(ZstdDecoder::new(file)?)
            }
            _ => {
                Box::new(BufReader::new(file))
            }
        };
        let reader = BufReader::new(reader);
        let mut ids = Vec::new();
        let mut data = Vec::with_capacity(100 * 1024 * 1024); 
        let mut offsets = vec![0];

        let mut current_id = String::new();
        let mut in_record = false;
        for line_res in reader.lines() {
            let line = line_res.context("Fail to read line")?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with(">") {
                if in_record {
                    ids.push(current_id);
                    offsets.push(data.len());
                }
            
                current_id = line[1..].split_whitespace().next().unwrap_or("unknown").to_string();
                
                in_record = true;
        } else{
            data.extend_from_slice(line.as_bytes());
        }
    }
    if in_record {
        ids.push(current_id);
        offsets.push(data.len());
    }
    println!("Loaded {} proteins from {:?}", ids.len(), path);
        
        Ok(Database {
            ids,
            data,
            offsets,
        })
    }
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
    pub fn get(&self, index: usize) -> Option<(&str, &[u8])> {
        if index < self.ids.len() {
            let start = self.offsets[index];
            let end = self.offsets[index + 1];
            Some((&self.ids[index], &self.data[start..end]))
        } else {
            None
        }
    }
}


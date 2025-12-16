use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::thread::current;
use anyhow::{Context, Result};

pub struct Database {
    pub ids: Vec<String>,
    pub sequences: Vec<Vec<u8>>,
}

impl Database {
    pub fn load_from_fasta<P: AsRef<Path>>(path:P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open file: {:?}", path))?;
        let reader = BufReader::new(file);

        let mut ids = Vec::new();
        let mut sequences = Vec::new();

        let mut current_id = String::new();
        let mut current_sequence = Vec::new();
        let mut in_record = false;
        for line_res in reader.lines() {
            let line = line_res.context("Fail to read line")?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with(">") {
                if in_record {
                    if current_id.is_empty() {
                        current_id= format!("unknown_{}", ids.len());
                    }
                    ids.push(current_id);
                    sequences.push(current_sequence);
                }
                current_id = line[1..].to_string();
                current_sequence = Vec::new();
                in_record = true;
        } else{
            current_sequence.extend_from_slice(line.as_bytes());
        }
    }
    if in_record && !current_sequence.is_empty() {
        ids.push(current_id);
        sequences.push(current_sequence);
    }
    println!("Loaded {} proteins from {:?}", ids.len(), path);
        
        Ok(Database {
            ids,
            sequences,
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
            Some((&self.ids[index], &self.sequences[index]))
        } else {
            None
        }
    }
}


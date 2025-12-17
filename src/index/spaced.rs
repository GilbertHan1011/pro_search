use rustc_hash::FxHashMap;
use crate::core::database::Database;
use crate::core::alphabet::{encode_spaced};
use crate::index::kmer::{ProteinId, Position};


pub struct SpacedIndex {
    pub map: FxHashMap<u64, Vec<(ProteinId, Position)>>,
    pub pattern: String,
    pub mask: Vec<bool>,
    pub weight: usize,
}

impl SpacedIndex {
    pub fn build(db: &Database, pattern: &str) -> Self {
        let mask: Vec<bool> = pattern.chars().map(|c| c == '1').collect();
        let weight = mask.iter().filter(|&&x| x).count();
        let span = mask.len();
        let mut map: FxHashMap<u64, Vec<(ProteinId, Position)>> = FxHashMap::default();

        println!("Building Spaced Index (Pattern: {}, Weight: {}, Span: {})...", pattern, weight, span);
        for (pid, seq) in db.sequences.iter().enumerate() {
            if seq.len() < span {
                continue;
            }

            for i in 0..=(seq.len() - span) {
                if let Some(encoded) = encode_spaced(&seq[i..i+span], &mask) {
                    map.entry(encoded)
                        .or_default()
                        .push((pid as ProteinId, i as Position));
                }
            }
        }
        Self {
            map,
            pattern: pattern.to_string(),
            mask: mask,
            weight,
        }
    }

    pub fn search_basic(&self, query_seq: &[u8], top_n: usize) -> Vec<(ProteinId, u32)> {
        let span = self.mask.len();
        let mut scores: FxHashMap<ProteinId, u32> = FxHashMap::default();
        for window in query_seq.windows(span) {
            if let Some(encoded) = encode_spaced(window, &self.mask) {
                if let Some(hits) = self.map.get(&encoded) {
                    for &(pid, _pos) in hits {
                        *scores.entry(pid).or_insert(0) += 1;
                    }
                }
            }
        }
        let mut candidates: Vec<(ProteinId, u32)> = scores.into_iter().collect();
        candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        if candidates.len() > top_n {
            candidates.truncate(top_n);
        }
        candidates
    }
    
        
}

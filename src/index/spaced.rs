use rustc_hash::FxHashMap;
use crate::core::database::Database;
use crate::core::alphabet::{encode_spaced};
use crate::index::kmer::{ProteinId, Position};
use smallvec::SmallVec;

pub type PostingsList = SmallVec<[(ProteinId, Position); 2]>;

pub struct SpacedIndex {
    pub map: FxHashMap<u64, PostingsList>,
    pub pattern: String,
    pub mask: Vec<bool>,
    pub weight: usize,
    pub num_proteins: usize,
}

impl SpacedIndex {
    pub fn build(db: &Database, pattern: &str) -> Self {
        let mask: Vec<bool> = pattern.chars().map(|c| c == '1').collect();
        let weight = mask.iter().filter(|&&x| x).count();
        let span = mask.len();
        let estimated_capacity = db.len() * 100;
        let mut map: FxHashMap<u64, PostingsList> = FxHashMap::with_capacity_and_hasher(estimated_capacity, Default::default());

        println!("Building Spaced Index (Pattern: {}, Weight: {}, Span: {})...", pattern, weight, span);
        for i in 0..db.len() {
            let pid = i as ProteinId;
            let start = db.offsets[i];
            let end = db.offsets[i+1];
            let seq = &db.data[start..end];
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
            mask,
            weight,
            num_proteins: db.len(),
        }
    }

    pub fn search_basic(&self, query_seq: &[u8], top_n: usize) -> Vec<(ProteinId, u32)> {
        let span = self.mask.len();
        let mut scores: Vec<u32> = vec![0; self.num_proteins];
        let mut active_pids: Vec<ProteinId> = Vec::new();
        for window in query_seq.windows(span) {
            if let Some(encoded) = encode_spaced(window, &self.mask) {
                if let Some(hits) = self.map.get(&encoded) {
                    // SmallVec 遍历极其高效
                    for &(pid, _pos) in hits {
                        let idx = pid as usize;
                        if scores[idx] == 0 {
                            active_pids.push(pid);
                        }
                        scores[idx] += 1;
                    }
                }
            }
        }
        let mut candidates: Vec<(ProteinId, u32)> = active_pids
            .into_iter()
            .map(|pid| (pid, scores[pid as usize]))
            .collect();
        candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        if candidates.len() > top_n {
            candidates.truncate(top_n);
        }
        candidates
    }
}

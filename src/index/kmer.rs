use rustc_hash::FxHashMap;
use crate::core::database::Database;
use crate::core::alphabet::encode_kmer;
use std::mem;

pub type ProteinId = u32;
pub type Position = u16;

pub struct KmerIndex {
    pub map: FxHashMap<u64, Vec<(ProteinId, Position)>>, // pub for test purpose
    pub k: usize,
}

impl KmerIndex {
    pub fn new(k: usize) -> Self {
        Self {
            map: FxHashMap::default(),
            k,
        }
    }
    pub fn build(db: &Database, k: usize) -> Self {
        let total_len: usize = db.sequences.iter().map(|s| s.len()).sum();
        let mut map: FxHashMap<u64, Vec<(ProteinId, Position)>> = FxHashMap::with_capacity_and_hasher(total_len, Default::default());
        println!("Building index with k={} for {} proteins...", k, db.len());
        for (prot_id, seq) in db.sequences.iter().enumerate() {
            let pid = prot_id as ProteinId;
            if seq.len() < k {
                continue;
            }
            for (pos, window) in seq.windows(k).enumerate() {
                if let Some(encoded) = encode_kmer(window) {
                    map.entry(encoded)
                    .or_default()
                    .push((pid, pos as Position));
                }
            }
        }
        println!("Index built! Total unique k-mers: {}", map.len());
        KmerIndex { map, k }

    }
    pub fn query(&self, encoded_kmer: u64) -> Option<&Vec<(ProteinId, Position)>> {
        self.map.get(&encoded_kmer)
    }

    pub fn search_basic(&self, query_seq: &[u8], top_n: usize) -> Vec<(ProteinId, u32)> {
        let mut scores: FxHashMap<ProteinId, u32> = FxHashMap::default();
        if query_seq.len() >= self.k {
            for window in query_seq.windows(self.k) {
                if let Some(encoded) = encode_kmer(window) {
                    if let Some(hits) = self.map.get(&encoded) {
                        for &(pid, _pos) in hits {
                            *scores.entry(pid).or_insert(0) += 1;
                        }
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
    pub fn memory_usage(&self) -> usize {
        let mut total_bytes = 0;

        // 1. Stack
        total_bytes += mem::size_of::<Self>();
        // 2. HashMap skeleton size
        let map_cap = self.map.capacity();
        // Key size
        total_bytes += map_cap * mem::size_of::<u64>();
        // Value size
        total_bytes += map_cap * mem::size_of::<Vec<(ProteinId, Position)>>();
        // Control Bytes
        total_bytes += map_cap * 1; 
        // 3. Payload Deep Size
        let item_size = mem::size_of::<(ProteinId, Position)>();
        for postings in self.map.values() {
            total_bytes += postings.capacity() * item_size;
        }

        total_bytes
    }
    
}
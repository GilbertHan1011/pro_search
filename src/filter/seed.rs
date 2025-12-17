use rustc_hash::FxHashMap;
use crate::index::kmer::KmerIndex;
use crate::index::kmer::{ProteinId};
use crate::core::alphabet::encode_kmer;

#[derive(Debug, Clone)]
pub struct Candidate {
    pub id: ProteinId,
    pub score: usize,
    pub best_diagonal: i32,
}

pub fn find_candidate(
    kmer_index: &KmerIndex, 
    query_seq: &[u8],
    min_diagonal: usize
) -> Vec<Candidate>{
    let mut protein_hit : FxHashMap<ProteinId,Vec<i32>> = FxHashMap::default();

    // 1. Sliding windows and query
    if query_seq.len() > kmer_index.k {
        for (q_pos, window) in query_seq.windows(kmer_index.k).enumerate() {
            if let Some(encoded) = encode_kmer(window){
                if let Some(hits) = kmer_index.query(encoded){
                    for &(pid, t_pos) in hits{
                        let diagonal = t_pos as i32 - q_pos as i32;
                        protein_hit.entry(pid).or_default().push(diagonal);
                    }
                }
            }
        }
    }


    // 2. Calculate scores and best diagonals
    let mut candidate = Vec::new();
    for (pid, mut diagonals) in protein_hit {
        diagonals.sort_unstable();
        let mut max_diagonal = 0;
        let mut max_hit = 0;
        let mut diagonal_count = 0;
        let mut current_diagonal = diagonals[0];
        for &diagonal in diagonals.iter() {
            if diagonal == current_diagonal {
                diagonal_count += 1;
            } else{
                if diagonal_count > max_hit {
                    max_hit = diagonal_count;
                    max_diagonal = current_diagonal;
                }
                current_diagonal = diagonal;
                diagonal_count = 1;
            }
        }
        if diagonal_count > max_hit {
            max_hit = diagonal_count;
            max_diagonal = current_diagonal;
        }
        if max_hit >= min_diagonal {
            candidate.push(Candidate { id: pid, score: max_hit, best_diagonal: max_diagonal });
        }
    }
    candidate.sort_unstable_by(|a, b| b.score.cmp(&a.score));

    candidate

}

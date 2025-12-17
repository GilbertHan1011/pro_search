use crate::index::kmer::ProteinId;

pub struct BenchmarkResult {
    pub recall_at_1: f64,
    pub recall_at_10: f64,
    pub mrr: f64, // Mean Reciprocal Rank
    pub avg_candidates: f64,
    pub avg_time_ms: f64,
}

pub fn calculate_metrics(
    candidates_list: &[Vec<(ProteinId, u32)>], 
    truths: &[ProteinId],
    total_time_ms: f64
) -> BenchmarkResult {
    let mut hits_at_1 = 0;
    let mut hits_at_10 = 0;
    let mut reciprocal_ranks = 0.0;
    let mut total_candidates = 0;

    for (candidates, &truth_id) in candidates_list.iter().zip(truths.iter()) {
        total_candidates += candidates.len();

        if let Some(rank) = candidates.iter().position(|(pid, _)| *pid == truth_id) {
            if rank == 0 { hits_at_1 += 1; }
            if rank < 10 { hits_at_10 += 1; }
            
            // MRR: 1 / (rank + 1)
            reciprocal_ranks += 1.0 / (rank as f64 + 1.0);
        }
    }

    let n = truths.len() as f64;
    BenchmarkResult {
        recall_at_1: hits_at_1 as f64 / n,
        recall_at_10: hits_at_10 as f64 / n,
        mrr: reciprocal_ranks / n,
        avg_candidates: total_candidates as f64 / n,
        avg_time_ms: total_time_ms / n,
    }
}
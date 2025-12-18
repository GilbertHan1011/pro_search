use crate::core::database::Database;
use crate::index::kmer::ProteinId;
use crate::filter::seed::Candidate;
use crate::align::simd;

#[derive(Debug, Clone)]
pub struct ExtensionResult {
    pub score: i32,
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
}

pub struct Scoring {
    pub match_score: i32,
    pub mismatch_score: i32,
}

// Use simple score here, maybe BLOSUM62 in the future
impl Scoring {
    pub fn default() -> Self {
        Self{match_score: 1, mismatch_score: -1 }
    }
}

pub fn extend_ungapped(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
    x_drop: i32,
) -> ExtensionResult {
    let (right_score, right_q_end, right_t_end) = simd::extend_direction_simd(
        query, target, scoring, 
        q_start, t_start, 1, x_drop
    );
    let (left_score, left_q_start, left_t_start) = if q_start > 0 && t_start > 0 {
        simd::extend_direction_simd(
            query, target, scoring, 
            q_start - 1, t_start - 1, -1, x_drop
        )
    } else {
        (0, q_start, t_start)
    };
    ExtensionResult { 
        score:left_score + right_score, 
        q_start:left_q_start, 
        q_end:right_q_end,  
        t_start:left_t_start, 
        t_end:right_t_end }
}

// Wrap ungapped function
pub fn refine_ungapped(
    query: &[u8],
    candidates: &[Candidate],
    db: &Database,
    scoring: &Scoring,
    x_drop: i32,
    top_n: usize,
) -> Vec<(ProteinId, ExtensionResult)> {
    
    let mut hits = Vec::with_capacity(candidates.len());

    for cand in candidates {
        // Get Target sequence
        let (_, target_seq) = match db.get(cand.id as usize) {
            Some(entry) => entry,
            None => continue,
        };
        let diag = cand.best_diagonal;
        
        // If diag > 0 (T > Q): Q=0, T=diag
        // If diag < 0 (Q > T): Q=-diag, T=0
        let (q_start, t_start) = if diag >= 0 {
            (0, diag as usize)
        } else {
            ((-diag) as usize, 0)
        };

        // Safety check: prevent out of bounds
        if t_start >= target_seq.len() || q_start >= query.len() {
            continue;
        }

        // Extend
        let result = extend_ungapped(
            query, target_seq, scoring, 
            q_start, t_start, x_drop
        );

        hits.push((cand.id, result));
    }

    // Sort by Score Descending
    // Use unstable sort for faster
    hits.sort_unstable_by(|a, b| b.1.score.cmp(&a.1.score));

    // Truncate (Top N)
    if hits.len() > top_n {
        hits.truncate(top_n);
    }

    hits
}
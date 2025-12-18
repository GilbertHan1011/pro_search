use bio::alignment::pairwise::{self, Scoring};
use bio::alignment::Alignment;

pub fn align_sw(
    query: &[u8], 
    target: &[u8], 
    gap_open: i32, 
    gap_extend: i32,
    match_score: i32,
    mismatch_score: i32
) -> Alignment {
    let scoring = Scoring::from_scores(
        gap_open, gap_extend, 
        match_score, mismatch_score
    );

    let mut aligner = pairwise::Aligner::with_scoring(scoring);

    aligner.local(query, target)
}

/// Extract a window around a center position, returning a slice reference and offset.
/// This avoids unnecessary allocations by returning a slice instead of a Vec.
pub fn extract_window(
    full_seq: &[u8], 
    center_pos: usize, 
    radius: usize
) -> (&[u8], usize) {
    let start = center_pos.saturating_sub(radius);
    let end = std::cmp::min(full_seq.len(), center_pos + radius);
    (&full_seq[start..end], start)
}
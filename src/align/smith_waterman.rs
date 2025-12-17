use bio::alignment::pairwise::{self, Scoring};
use bio::alignment::Alignment;

pub fn align_sw(
    query: &[u8], 
    target: &[u8], 
    gap_open: i32, 
    gap_extend: i32
) -> Alignment {
    let scoring = Scoring::from_scores(5, -4, gap_open, gap_extend);

    let mut aligner = pairwise::Aligner::with_scoring(scoring);

    aligner.local(query, target)
}

pub fn extract_window(
    full_seq: &[u8], 
    center_pos: usize, 
    radius: usize
) -> (Vec<u8>, usize) {
    let start = center_pos.saturating_sub(radius);

    let end = std::cmp::min(full_seq.len(), center_pos + radius);

    (full_seq[start..end].to_vec(), start)
}
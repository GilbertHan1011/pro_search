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

fn extend_direction(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
    step: isize,
    x_drop: i32,
) -> (i32, usize, usize) {
    let mut best_score = 0;
    let mut current_score = 0;
    let mut best_q_idx = q_start;
    let mut best_t_idx = t_start;

    let mut q_idx = q_start as isize;
    let mut t_idx = t_start as isize;

    loop {
        if q_idx < 0 || t_idx < 0 || q_idx as usize >= query.len()|| t_idx as usize >= target.len() {
            break;
        }

        let q_char = query[q_idx as usize];
        let t_char = target[t_idx as usize];

        if q_char == t_char {
            current_score += scoring.match_score;
        } else {
            current_score += scoring.mismatch_score;
        }
        
        if current_score > best_score {
            best_score = current_score;
            best_q_idx = q_idx as usize;
            best_t_idx = t_idx as usize;
        } else if current_score < best_score - x_drop {
            break;
        }

        q_idx += step;
        t_idx += step;
    }

    (best_score, best_q_idx, best_t_idx)
}

pub fn extend_ungapped(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
    x_drop: i32,
) -> ExtensionResult {
    let (right_score, right_q_end, right_t_end) = extend_direction(query, target, 
        scoring, q_start, 
        t_start, 1, 
        x_drop);
    let (left_score, left_q_start, left_t_start) = if q_start > 0 && t_start > 0 {
        extend_direction(query, target, 
        scoring, q_start -1,  // avoid double counting
        t_start -1, -1, 
        x_drop)
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
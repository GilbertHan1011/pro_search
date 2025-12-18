#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;

use crate::align::ungapped::Scoring;

/// Optimized right-extension using AVX2.
/// Falls back to scalar if AVX2 is not available or for left-extension.
pub fn extend_direction_simd(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
    step: isize,
    x_drop: i32,
) -> (i32, usize, usize) {
    // Only optimize Right Extension (+1) with AVX2
    if step == 1 && is_x86_feature_detected!("avx2") {
        unsafe { extend_right_avx2(query, target, scoring, q_start, t_start, x_drop) }
    } else {
        // Fallback for Left Extension (-1) or non-AVX CPUs
        extend_direction_scalar(query, target, scoring, q_start, t_start, step, x_drop)
    }
}

/// The Scalar Fallback (Your original function)
fn extend_direction_scalar(
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
        // Boundary check
        if q_idx < 0 || t_idx < 0 || q_idx as usize >= query.len() || t_idx as usize >= target.len() {
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

/// AVX2 Optimized Implementation for Step = 1
#[target_feature(enable = "avx2")]
unsafe fn extend_right_avx2(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
    x_drop: i32,
) -> (i32, usize, usize) {
    let mut best_score = 0;
    let mut current_score = 0;
    let mut best_q_idx = q_start;
    let mut best_t_idx = t_start;

    let mut q_curr = q_start;
    let mut t_curr = t_start;
    
    let q_len = query.len();
    let t_len = target.len();

    // 1. SIMD Loop: Process 32 bytes at a time
    while q_curr + 32 <= q_len && t_curr + 32 <= t_len {
        // Load 32 bytes from query and target (Unaligned load is fast on modern CPUs)
        let v_q = _mm256_loadu_si256(query.as_ptr().add(q_curr) as *const __m256i);
        let v_t = _mm256_loadu_si256(target.as_ptr().add(t_curr) as *const __m256i);

        // Compare: result is 0xFF for equal, 0x00 for not equal
        let v_cmp = _mm256_cmpeq_epi8(v_q, v_t);

        // Create Mask: Collapses 32 bytes into a 32-bit integer
        // Bit i is 1 if match, 0 if mismatch
        let mask = _mm256_movemask_epi8(v_cmp) as u32;

        // Optimization: Fast Path
        // Calculate the score change for the entire block
        let matches = mask.count_ones() as i32;
        let mismatches = 32 - matches;
        let block_delta = matches * scoring.match_score + mismatches * scoring.mismatch_score;

        // "Lookahead" Check:
        // If adding the whole block doesn't trigger X-drop AND ensures a new best score,
        // we can skip the bit-by-bit check.
        // (Note: This is a heuristic. For strict X-drop, we usually iterate bits. 
        //  Here we iterate bits to be safe and accurate).
        
        // --- Bit Iteration (No Memory Access) ---
        // This loop runs entirely in registers, very fast.
        let mut temp_mask = mask;
        for i in 0..32 {
            // Check lowest bit
            if (temp_mask & 1) != 0 {
                current_score += scoring.match_score;
            } else {
                current_score += scoring.mismatch_score;
            }

            if current_score > best_score {
                best_score = current_score;
                best_q_idx = q_curr + i;
                best_t_idx = t_curr + i;
            } else if current_score < best_score - x_drop {
                // Drop detected inside the block!
                return (best_score, best_q_idx, best_t_idx);
            }
            
            // Shift to next bit
            temp_mask >>= 1;
        }

        q_curr += 32;
        t_curr += 32;
    }

    // 2. Tail Loop: Process remaining characters (< 32)
    // We pass the current state to the scalar function to finish up
    if q_curr < q_len && t_curr < t_len {
        // We perform a small "mini-scalar" extension here manually 
        // to avoid reconstructing the whole loop state wrapper
        let mut q_i = q_curr;
        let mut t_i = t_curr;
        while q_i < q_len && t_i < t_len {
            let val_q = *query.get_unchecked(q_i);
            let val_t = *target.get_unchecked(t_i);

            if val_q == val_t {
                current_score += scoring.match_score;
            } else {
                current_score += scoring.mismatch_score;
            }

            if current_score > best_score {
                best_score = current_score;
                best_q_idx = q_i;
                best_t_idx = t_i;
            } else if current_score < best_score - x_drop {
                break;
            }
            q_i += 1;
            t_i += 1;
        }
    }

    (best_score, best_q_idx, best_t_idx)
}
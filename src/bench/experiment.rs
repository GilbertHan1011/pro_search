use std::time::Instant;
use crate::core::database::Database;
use crate::index::kmer::{KmerIndex, ProteinId};
use crate::align::ungapped::{self, Scoring};
use crate::bench::query_gen::{self, QueryConfig};
use crate::bench::metric::{self, BenchmarkResult};

use crate::filter::seed; // Step 2
use crate::index::spaced::SpacedIndex; // Step 6
use crate::align::smith_waterman; // Step 5

struct ExpResult {
    name: String,
    recall_10: f64,
    mrr: f64,
    avg_time: f64,
    candidates: f64,
}

impl ExpResult {
    fn print(&self) {
        println!("{:20} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms | Cands: {:.1}", 
                self.name, self.recall_10, self.mrr, self.avg_time, self.candidates);
    }
}

fn print_comparison(title: &str, baseline: &BenchmarkResult, refined: &BenchmarkResult, label_base: &str, label_refined: &str) {
    println!("\n=== {} ===", title);
    println!("{:<20} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms", label_base, baseline.recall_at_10, baseline.mrr, baseline.avg_time_ms);
    println!("{:<20} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms", label_refined, refined.recall_at_10, refined.mrr, refined.avg_time_ms);
    
    let mrr_gain = (refined.mrr - baseline.mrr) / baseline.mrr * 100.0;
    println!(">> MRR Improvement: {:.2}%", mrr_gain);
}

pub fn run_k_tradeoff(db: &Database) {
    println!("\n=== Task 1: K-mer Trade-off Analysis ===");

    let config = QueryConfig { length: 60, sub_rate: 0.10, indel_rate: 0.0 };
    let queries = query_gen::sample_queries(db, 100, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    for k in 3..=7 {
        let start = Instant::now();
        let index = KmerIndex::build(db, k);
        let build_time = start.elapsed().as_millis();


        let start_search = Instant::now();
        let mut candidates_list = Vec::new();
        
        for q in &queries {
            let hits = index.search_basic(&q.sequence, 10);
            candidates_list.push(hits);
        }
        let total_time = start_search.elapsed().as_millis() as f64;

        let metrics = metric::calculate_metrics(&candidates_list, &truths, total_time);
        
        let res = ExpResult {
            name: format!("k={}", k),
            recall_10: metrics.recall_at_10,
            mrr: metrics.mrr,
            avg_time: metrics.avg_time_ms,
            candidates: metrics.avg_candidates,
        };
        res.print();
        println!("   (Index Build Time: {} ms, Size: {} k-mers)", build_time, index.map.len());
    }
}

pub fn run_filter_comparison(db: &Database) {
    println!("\n=== Task 2: Diagonal Filtering vs Voting (k=5) ===");
    
    let config = QueryConfig { length: 60, sub_rate: 0.20, indel_rate: 0.0 };
    let queries = query_gen::sample_queries(db, 100, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    let k = 5;
    let index = KmerIndex::build(db, k);

    let start_a = Instant::now();
    let res_a: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index.search_basic(&q.sequence, 10))
        .collect();
    let metrics_a = metric::calculate_metrics(&res_a, &truths, start_a.elapsed().as_millis() as f64);
    
    ExpResult {
        name: "Method A (Voting)".to_string(),
        recall_10: metrics_a.recall_at_10,
        mrr: metrics_a.mrr,
        avg_time: metrics_a.avg_time_ms,
        candidates: metrics_a.avg_candidates
    }.print();

    let start_b = Instant::now();
    let res_b: Vec<Vec<(ProteinId, u32)>> = queries.iter().map(|q| {
        let cands = seed::find_candidate(&index, &q.sequence, 2);
        cands.into_iter().take(10).map(|c| (c.id, c.score as u32)).collect()
    }).collect();
    let metrics_b = metric::calculate_metrics(&res_b, &truths, start_b.elapsed().as_millis() as f64);

    ExpResult {
        name: "Method B (Diag)".to_string(),
        recall_10: metrics_b.recall_at_10,
        mrr: metrics_b.mrr,
        avg_time: metrics_b.avg_time_ms,
        candidates: metrics_b.avg_candidates
    }.print();
}

pub fn run_spaced_seed_test(db: &Database) {
    println!("\n=== Task 4: Spaced Seeds vs Contiguous (High Mutation) ===");
    
    // 极高难度：30% 突变率
    let config = QueryConfig { length: 60, sub_rate: 0.30, indel_rate: 0.0 };
    let queries = query_gen::sample_queries(db, 200, &config); // 跑多一点样本
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    // 1. Contiguous k=5 (Weight=5, Span=5)
    let index_k5 = KmerIndex::build(db, 5);
    
    // 2. Spaced Pattern (Weight=5, Span=7)
    // 1101011 -> Weight 5
    let index_spaced = SpacedIndex::build(db, "1101011");

    // --- Run Contiguous ---
    let start_1 = Instant::now();
    let res_1: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index_k5.search_basic(&q.sequence, 10))
        .collect();
    let m1 = metric::calculate_metrics(&res_1, &truths, start_1.elapsed().as_millis() as f64);
    
    ExpResult {
        name: "Contiguous (11111)".to_string(),
        recall_10: m1.recall_at_10,
        mrr: m1.mrr,
        avg_time: m1.avg_time_ms,
        candidates: m1.avg_candidates
    }.print();

    // --- Run Spaced ---
    let start_2 = Instant::now();
    let res_2: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index_spaced.search_basic(&q.sequence, 10))
        .collect();
    let m2 = metric::calculate_metrics(&res_2, &truths, start_2.elapsed().as_millis() as f64);

    ExpResult {
        name: "Spaced (1101011)".to_string(),
        recall_10: m2.recall_at_10,
        mrr: m2.mrr,
        avg_time: m2.avg_time_ms,
        candidates: m2.avg_candidates
    }.print();
}


/// Task 3 / Step 5 Test: Indel Robustness & SW Refinement
pub fn run_indel_test(db: &Database) {
    println!("\nGenerating Indel Queries...");

    let config = QueryConfig { length: 80, sub_rate: 0.05, indel_rate: 0.10 };
    let queries = query_gen::sample_queries(db, 50, &config); // 跑 50 个样本
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    let k = 5;
    let index = KmerIndex::build(db, k);
    let scoring = Scoring::default();

    // 存储两种策略的结果
    let mut res_ungapped = Vec::new(); // Baseline
    let mut res_sw = Vec::new();       // Refined

    let start_total = Instant::now();

    for q in &queries {
        // --- Step 2: Seeding ---
        let candidates = seed::find_candidate(&index, &q.sequence, 2);
        
        // --- Step 3: Ungapped Extension (Baseline) ---
        let mut ungapped_hits = Vec::new();
        for cand in candidates.iter().take(50) {
            let (_, target_seq) = db.get(cand.id as usize).unwrap();

            let q_start = 0; 
            let t_start = (cand.best_diagonal + 0).max(0) as usize; 

            let ext = ungapped::extend_ungapped(
                &q.sequence, &target_seq, &scoring,
                q_start, t_start, k, 10
            );
            ungapped_hits.push((cand.id, ext));
        }
        
        ungapped_hits.sort_by(|a, b| b.1.score.cmp(&a.1.score));
        
        res_ungapped.push(ungapped_hits.iter().map(|(id, ext)| (*id, ext.score as u32)).collect());

        // --- Step 5: Smith-Waterman Refinement (Advanced) ---
        let mut sw_hits = Vec::new();
        let window_radius = 50;

        for (id, ext) in ungapped_hits.iter().take(20) {
            let (_, target_full) = db.get(*id as usize).unwrap();

            let q_center = (ext.q_start + ext.q_end) / 2;
            let t_center = (ext.t_start + ext.t_end) / 2;
            
            // extract_window
            let (q_sub, _) = smith_waterman::extract_window(&q.sequence, q_center, window_radius);
            let (t_sub, _) = smith_waterman::extract_window(target_full, t_center, window_radius);

            // 2. Smith-Waterman
            // Gap Open = -10, Gap Extend = -1
            let align = smith_waterman::align_sw(&q_sub, &t_sub, -10, -1);
            
            sw_hits.push((*id, align.score));
        }

        sw_hits.sort_by(|a, b| b.1.cmp(&a.1));
        res_sw.push(sw_hits.into_iter().map(|(id, s)| (id, s as u32)).collect());
    }

    let time_per_query = start_total.elapsed().as_millis() as f64 / 50.0;

    let m_ungapped = metric::calculate_metrics(&res_ungapped, &truths, time_per_query);
    let m_sw = metric::calculate_metrics(&res_sw, &truths, time_per_query); // 注意：这里时间其实是累加的

    print_comparison("Step 5 Test: Indel Robustness (10% Indels)", &m_ungapped, &m_sw, "Ungapped Only", "Ungapped + SW");
}



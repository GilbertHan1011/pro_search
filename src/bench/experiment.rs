use std::time::Instant;
use crate::core::database::Database;
use crate::index::kmer::{KmerIndex, ProteinId};
use crate::align::ungapped::{Scoring, refine_ungapped};
use crate::bench::query_gen::{self, QueryConfig};
use crate::bench::helper::{print_comparison, run_gapped_wrapper,create_csv_writer,write_metric_csv};
use crate::bench::metric::calculate_metrics;
use crate::filter::seed; // Step 2
use crate::index::spaced::SpacedIndex; // Step 6
use crate::align::smith_waterman; // Step 5
use std::path::Path;
use std::io::Write;
pub struct ExpResult {
    name: String,
    recall_1: f64,
    recall_10: f64,
    mrr: f64,
    avg_time: f64,
    candidates: f64,
}

impl ExpResult {
    fn print(&self) {
        println!("{:20} | R@1: {:.4} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms | Cands: {:.1}", 
                self.name, self.recall_1, self.recall_10, self.mrr, self.avg_time, self.candidates);
    }
    fn csv_header() -> &'static str {
        "Name,Recall_1,Recall_10,MRR,Avg_Time_ms,Avg_Candidates"
    }
    fn to_csv_line(&self) -> String {
        format!("{},{:.6},{:.6},{:.6},{:.6},{:.2}\n", 
                self.name, self.recall_1, self.recall_10, self.mrr, self.avg_time, self.candidates)
    }
}


pub fn run_k_tradeoff(
    db: &Database,top_n: usize, mutate: bool, 
    length: usize, sub_rate: f64, 
    indel_rate: f64,sample_num: usize,
    csv_path: Option<&Path>
) {
    println!("\n=== Task 1: K-mer Trade-off Analysis ===");
    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "{}", ExpResult::csv_header()).unwrap();
    }
    let config = if mutate {
        QueryConfig { length, sub_rate, indel_rate }
    } else {
        QueryConfig { length, sub_rate: 0.0, indel_rate: 0.0 }
    };
    let queries = query_gen::sample_queries(db, sample_num, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    for k in 3..=7 {
        let start = Instant::now();
        let index = KmerIndex::build(db, k);
        let build_time = start.elapsed().as_millis();


        let start_search = Instant::now();
        let mut candidates_list = Vec::new();
        
        for q in &queries {
            let hits = index.search_basic(&q.sequence, top_n);
            candidates_list.push(hits);
        }
        let total_time = start_search.elapsed().as_millis() as f64;

        let metrics = calculate_metrics(&candidates_list, &truths, total_time);
        
        let res = ExpResult {
            name: format!("k={}", k),
            recall_1: metrics.recall_at_1,
            recall_10: metrics.recall_at_10,
            mrr: metrics.mrr,
            avg_time: metrics.avg_time_ms,
            candidates: metrics.avg_candidates,
        };
        res.print();
        println!("   (Index Build Time: {} ms, Size: {} k-mers)", build_time, index.map.len());
        if let Some(w) = &mut csv_writer {
            write!(w, "{}", res.to_csv_line()).unwrap();
        }
    }
}

pub fn run_filter_comparison(
    db: &Database,top_n: usize, k: usize, mutate: bool, 
    length: usize, sub_rate: f64, indel_rate: f64,sample_num: usize,
    csv_path: Option<&Path>
) {
    println!("\n=== Task 2: Diagonal Filtering vs Voting (k={}) ===", k);
    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "{}", ExpResult::csv_header()).unwrap();
    }
    let config = if mutate {
        QueryConfig { length, sub_rate, indel_rate }
    } else {
        QueryConfig { length, sub_rate: 0.0, indel_rate: 0.0 }
    };
    let queries = query_gen::sample_queries(db, sample_num, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    let index = KmerIndex::build(db, k);

    let start_a = Instant::now();
    let res_a: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index.search_basic(&q.sequence, top_n))
        .collect();
    let metrics_a = calculate_metrics(&res_a, &truths, start_a.elapsed().as_millis() as f64);
    
    let result_a = ExpResult {
        name: "Method A (Voting)".to_string(),
        recall_1: metrics_a.recall_at_1,
        recall_10: metrics_a.recall_at_10,
        mrr: metrics_a.mrr,
        avg_time: metrics_a.avg_time_ms,
        candidates: metrics_a.avg_candidates
    };
    result_a.print();
    if let Some(w) = &mut csv_writer {
        write!(w, "{}", result_a.to_csv_line()).unwrap();
    }

    let start_b = Instant::now();
    let res_b: Vec<Vec<(ProteinId, u32)>> = queries.iter().map(|q| {
        let cands = seed::find_candidate(&index, &q.sequence, 2);
        cands.into_iter().take(10).map(|c| (c.id, c.score as u32)).collect()
    }).collect();
    let metrics_b = calculate_metrics(&res_b, &truths, start_b.elapsed().as_millis() as f64);

    let result_b = ExpResult {
        name: "Method B (Diag)".to_string(),
        recall_1: metrics_b.recall_at_1,
        recall_10: metrics_b.recall_at_10,
        mrr: metrics_b.mrr,
        avg_time: metrics_b.avg_time_ms,
        candidates: metrics_b.avg_candidates
    };
    result_b.print();
    if let Some(w) = &mut csv_writer {
        write!(w, "{}", result_b.to_csv_line()).unwrap();
    }
}

pub fn run_ungapped_test(
    db: &Database, sample_num: usize, 
    top_n: usize, k: usize, 
    length: usize, sub_rate: f64, 
    indel_rate: f64,x_drop: usize,
    csv_path: Option<&Path>
) {

    println!("\n=== Task 3: Ungapped Extension ===");
    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "Scenario,Recall_1,Recall_10,MRR,Avg_Time_ms,Avg_Candidates").unwrap();
    }
    let config_sub  = QueryConfig{length, sub_rate, indel_rate:0.0};
    let config_indel = QueryConfig{length, sub_rate:0.0, indel_rate};
    
    let index = KmerIndex::build(db, k);
    let scoring = Scoring::default();
    let metrics_sub = run_gapped_wrapper(db,&index,top_n,sample_num,&config_sub,&scoring,x_drop as i32);
    let metrics_indel = run_gapped_wrapper(db,&index,top_n,sample_num,&config_indel,&scoring,x_drop as i32);


    print_comparison("Ungapped Extension", &metrics_sub, 
    &metrics_indel, "Sub Only", "Indel only");
    if let Some(w) = &mut csv_writer {
        write_metric_csv(w, "Sub Only", &metrics_sub);
        write_metric_csv(w, "Indel Only", &metrics_indel);
    }
}


pub fn run_stress_all(db: &Database, sample_num: usize, top_n: usize, csv_path: Option<&Path>) {
    println!("\n===============================================================");
    println!("   STRESS TEST: K-mer Tradeoffs & Mutation Robustness");
    println!("===============================================================");
    println!("{:<4} | {:<8} | {:<8} | {:<8} | {:<8} | {:<8}| {:<8} | {:<10}", 
            "K", "Mut_Rate", "R@1", "R@10", "MRR", "Time(ms)", "Cands", "Mem(MB)");
    println!("---------------------------------------------------------------");

    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "K,Mutation_Rate,Recall_1,Recall_10,MRR,Avg_Time_ms,Avg_Candidates,Memory_MB").unwrap();
    }

    let k_values = vec![3, 4, 5, 6, 7, 8];
    let mutation_rates = vec![0.0, 0.10, 0.20, 0.30];

    for &sub_rate in &mutation_rates {
        let config = QueryConfig { 
            length: 60, 
            sub_rate, 
            indel_rate: 0.0 // Focus on Substitution
        };
        let queries = query_gen::sample_queries(db, sample_num, &config);
        let truths: Vec<_> = queries.iter().map(|q| q.original_pid).collect();

        for &k in &k_values {
            let index = KmerIndex::build(db, k);
            let mem_bytes = index.memory_usage();
            let mem_mb = mem_bytes as f64 / 1024.0 / 1024.0;

            let start_search = Instant::now();
            let mut candidates_list = Vec::new();

            for q in &queries {
                // Use Step 2 (Diagonal Filter) or Step 1 (Basic)
                let hits = index.search_basic(&q.sequence, top_n);
                candidates_list.push(hits);
            }
            
            let total_time = start_search.elapsed().as_millis() as f64;
            
            let m = calculate_metrics(&candidates_list, &truths, total_time);

            println!("{:<4} | {:<8.1} | {:<8.4} | {:<8.4} | {:<8.4} | {:<8.1} | {:<10.2} | {:<10.2}", 
                    k, sub_rate, m.recall_at_1, m.recall_at_10, m.mrr, m.avg_time_ms, m.avg_candidates, mem_mb);
            if let Some(w) = &mut csv_writer {
                write!(w, "{},{:.6},{:.6},{:.6},{:.6},{:.2},{:.2},{:.2}\n", 
                    k, sub_rate, m.recall_at_1, m.recall_at_10, m.mrr, m.avg_time_ms, m.avg_candidates, mem_mb).unwrap();
            }
    }
    println!("---------------------------------------------------------------");
}
}


/// Task 5 : Indel Robustness & SW Refinement
pub fn run_indel_test(
    db: &Database, sample_num: usize, 
    top_n: usize, k: usize, 
    length: usize, sub_rate: f64, 
    indel_rate: f64, x_drop: usize,
    csv_path: Option<&Path>
) {
    println!("\n=== Task 5: Indel Robustness & SW Refinement ===");
    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "Scenario,Recall_1,Recall_10,MRR,Avg_Time_ms,Avg_Candidates").unwrap();
    }
    let config = QueryConfig { length, sub_rate, indel_rate };
    let queries = query_gen::sample_queries(db, sample_num, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    let index = KmerIndex::build(db, k);
    let scoring = Scoring::default();

    let mut res_ungapped = Vec::new(); // Baseline
    let mut res_sw = Vec::new();       // Refined

    let start_total = Instant::now();

    for q in &queries {
        // --- Step 2: Seeding ---
        let candidates = seed::find_candidate(&index, &q.sequence, 2);
        // --- Step 3: Ungapped Extension (Baseline) ---
        let ungapped_hits = refine_ungapped(&q.sequence, &candidates, db, &scoring, x_drop as i32, top_n);

        res_ungapped.push(ungapped_hits.iter().map(|(id, ext)| (*id, ext.score as u32)).collect());

        // --- Step 5: Smith-Waterman Refinement (Advanced) ---
        let mut sw_hits = Vec::new();
        let window_radius = 80;

        for (id, ext) in ungapped_hits.iter() {
            let (_, target_full) = db.get(*id as usize).unwrap();

            let q_center = (ext.q_start + ext.q_end) / 2;
            let t_center = (ext.t_start + ext.t_end) / 2;
            
            // extract_window
            let (q_sub, _) = smith_waterman::extract_window(&q.sequence, q_center, window_radius);
            let (t_sub, _) = smith_waterman::extract_window(target_full, t_center, window_radius);

            // 2. Smith-Waterman
            // Gap Open = -10, Gap Extend = -1
            let align = smith_waterman::align_sw(&q_sub, &t_sub, -10, -1, 1, -1);
            
            sw_hits.push((*id, align.score));
        }

        sw_hits.sort_by(|a, b| b.1.cmp(&a.1));
        res_sw.push(sw_hits.into_iter().map(|(id, s)| (id, s as u32)).collect());
    }

    let time_per_query = start_total.elapsed().as_millis() as f64 / 50.0;

    let m_ungapped = calculate_metrics(&res_ungapped, &truths, time_per_query);
    let m_sw = calculate_metrics(&res_sw, &truths, time_per_query); 

    print_comparison("Step 5 Test: Indel Robustness (10% Indels)", &m_ungapped, 
    &m_sw, "Ungapped Only", "Ungapped + SW");
    if let Some(w) = &mut csv_writer {
        write_metric_csv(w, "Ungapped Only", &m_ungapped);
        write_metric_csv(w, "Ungapped + SW", &m_sw);
    }
}


pub fn run_spaced_seed_test(
    db: &Database, top_n: usize,
    k: usize,pattern: &str,sample_num: usize, 
    sub_rate: f64, indel_rate: f64,length: usize,
    csv_path: Option<&Path>
) {
    println!("\n=== Task 6: Spaced Seeds vs Contiguous (High Mutation) ===");
    let mut csv_writer = create_csv_writer(csv_path);
    if let Some(w) = &mut csv_writer {
        writeln!(w, "{}", ExpResult::csv_header()).unwrap();
    }
    // 30% mutation rate
    let config = QueryConfig { length, sub_rate, indel_rate };
    let queries = query_gen::sample_queries(db, sample_num, &config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    // 1. Contiguous k=5 (Weight=5, Span=5)
    let index_k5 = KmerIndex::build(db, k);
    
    // 2. Spaced Pattern (Weight=5, Span=7)
    // 1101011 -> Weight 5
    let index_spaced = SpacedIndex::build(db, pattern);

    if index_spaced.weight != k {
        panic!("Spaced index weight does not match k");
    }

    // --- Run Contiguous ---
    let start_1 = Instant::now();
    let res_1: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index_k5.search_basic(&q.sequence, top_n))
        .collect();
    let m1 = calculate_metrics(&res_1, &truths, start_1.elapsed().as_millis() as f64);
    
    let result_1 = ExpResult {
        name: "Contiguous (11111)".to_string(),
        recall_1: m1.recall_at_1,
        recall_10: m1.recall_at_10,
        mrr: m1.mrr,
        avg_time: m1.avg_time_ms,
        candidates: m1.avg_candidates
    };
    result_1.print();
    if let Some(w) = &mut csv_writer {
        write!(w, "{}", result_1.to_csv_line()).unwrap();
    }

    // --- Run Spaced ---
    let start_2 = Instant::now();
    let res_2: Vec<Vec<(ProteinId, u32)>> = queries.iter()
        .map(|q| index_spaced.search_basic(&q.sequence, top_n))
        .collect();
    let m2 = calculate_metrics(&res_2, &truths, start_2.elapsed().as_millis() as f64);

    let result_2 = ExpResult {
        name: "Spaced (1101011)".to_string(),
        recall_1: m2.recall_at_1,
        recall_10: m2.recall_at_10,
        mrr: m2.mrr,
        avg_time: m2.avg_time_ms,
        candidates: m2.avg_candidates
    };
    result_2.print();
    if let Some(w) = &mut csv_writer {
        write!(w, "{}", result_2.to_csv_line()).unwrap();
    }
}

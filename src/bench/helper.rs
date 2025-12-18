//Helper function like printing and wrapper
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use crate::bench::metric::{BenchmarkResult,calculate_metrics};
use crate::core::database::Database;
use crate::bench::query_gen::sample_queries;
use crate::filter::seed::find_candidate;
use crate::align::ungapped::{refine_ungapped,Scoring};
use crate::bench::query_gen::QueryConfig;
use std::time::Instant;
use crate::index::kmer::{ProteinId, KmerIndex};


pub fn print_comparison(title: &str, baseline: &BenchmarkResult, refined: &BenchmarkResult, label_base: &str, label_refined: &str) {
    println!("\n=== {} ===", title);
    println!("{:<20} | R@1: {:.4} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms", label_base, baseline.recall_at_1, baseline.recall_at_10, baseline.mrr, baseline.avg_time_ms);
    println!("{:<20} | R@1: {:.4} | R@10: {:.4} | MRR: {:.4} | Time: {:.4}ms", label_refined, refined.recall_at_1, refined.recall_at_10, refined.mrr, refined.avg_time_ms);
    
    let mrr_gain = (refined.mrr - baseline.mrr) / baseline.mrr * 100.0;
    println!(">> MRR Improvement: {:.2}%", mrr_gain);
}

pub fn create_csv_writer(path: Option<&Path>) -> Option<BufWriter<File>> {
    if let Some(p) = path {
        match File::create(p) {
            Ok(f) => Some(BufWriter::new(f)),
            Err(e) => {
                eprintln!("⚠️ Warning: Could not create CSV file {:?}: {}", p, e);
                None
            }
        }
    } else {
        None
    }
}

pub fn write_metric_csv<W: Write>(w: &mut W, name: &str, m: &crate::bench::metric::BenchmarkResult) {
    writeln!(w, "{},{:.6},{:.6},{:.6},{:.6},{:.2}", 
            name, m.recall_at_1, m.recall_at_10, m.mrr, m.avg_time_ms, m.avg_candidates).unwrap();
}

pub fn run_gapped_wrapper(
    db: &Database, index: &KmerIndex, top_n: usize, 
    sample_num: usize, config: &QueryConfig, 
    scoring: &Scoring, x_drop: i32) -> BenchmarkResult {
    let queries = sample_queries(db, sample_num, config);
    let truths: Vec<ProteinId> = queries.iter().map(|q| q.original_pid).collect();

    let mut results = Vec::with_capacity(queries.len());
    let start_time = Instant::now();

    for q in &queries {
        let candidates = find_candidate(index, &q.sequence, 2);

        // Step 3: Ungapped Extension 
        let refined_hits = refine_ungapped(
            &q.sequence, 
            &candidates, 
            db, 
            scoring, 
            x_drop, 
            top_n
        );

        let hits_formatted: Vec<(ProteinId, u32)> = refined_hits
            .into_iter()
            .map(|(id, ext)| (id, std::cmp::max(0, ext.score) as u32))
            .collect();
            
        results.push(hits_formatted);
    }

    let total_time = start_time.elapsed().as_millis() as f64;
    calculate_metrics(&results, &truths, total_time)
    
}

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use std::time::Instant;

use pro_search::core::database::Database;
use pro_search::index::kmer::KmerIndex;
use pro_search::index::spaced;
use pro_search::filter::{seed};
use pro_search::align::{ungapped, smith_waterman};
use pro_search::bench::experiment;


#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, value_name = "DB_FILE")]
    database: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Search {
        #[arg(long)]
        query: Option<String>,
        #[arg(long)]
        query_file: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = SearchMode::Auto)]
        mode: SearchMode,
        #[arg(short, long, default_value_t = 5)]
        k: usize,
        #[arg(short, long, default_value_t = 10)]
        n: usize,
    },
    Bench {
        #[arg(value_enum)] // Takes the enum as a required positional argument
        task: BenchTask,
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum SearchMode {
    /// K-mer Voting
    Basic,
    /// Diagonal Filtering
    Diagonal,
    /// Spaced Seeds
    Spaced,
    /// Full Pipeline (Index -> Filter -> Ungapped -> Smith-Waterman)
    Auto, 
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum BenchTask {
    /// Task 1: K-mer trade-off
    K,
    /// Task 2: Diagonal Filter vs Voting
    Filter,
    /// Task 3 & 5: Indel Robustness
    Indel,
    /// Task 4 & 6: Spaced Seeds
    Spaced,
    /// Run all benchmarks
    All,
}
fn main() {
    let args = Args::parse();
    if !args.database.exists() {
        eprintln!("❌ Error: Database file not found: {:?}", args.database);
        std::process::exit(1);
    }
    println!("Loading database from {:?}...", args.database);
    let start_load = Instant::now();
    let db = match Database::load_from_fasta(&args.database) {
        Ok(db) => db,
        Err(e) => {
            eprintln!("❌ Failed to load database: {}", e);
            std::process::exit(1);
        }
    };
    println!("✅ Database loaded in {:.2?} ({} proteins)", start_load.elapsed(), db.len());


    match args.command {
        Commands::Search { query, query_file, mode, k, n } => {
            // Collect all Queries
            let mut queries = Vec::new();
            
            // Source A: Command line input
            if let Some(q_str) = query {
                queries.push(("CommandLine_Query".to_string(), q_str.into_bytes()));
            }
            
            // Source B: File input
            if let Some(q_path) = query_file {
                match Database::load_from_fasta(&q_path) {
                    Ok(q_db) => {
                        for (id, seq) in q_db.ids.into_iter().zip(q_db.sequences.into_iter()) {
                            queries.push((id, seq));
                        }
                    },
                    Err(e) => eprintln!("⚠️ Warning: Failed to load query file: {}", e),
                }
            }
            if queries.is_empty() {
                eprintln!("❌ Error: No query provided. Use --query or --query-file.");
                return;
            }

            println!("Running search for {} queries (Mode: {:?}, k={})...", queries.len(), mode, k);
            
            // 根据模式构建索引 (只需构建一次)
            let start_idx = Instant::now();
            let index = KmerIndex::build(&db, k); // Basic/Diag/Auto 需要普通索引
            // 注意：如果是 Spaced 模式，我们需要另外构建 SpacedIndex
            let spaced_index = if mode == SearchMode::Spaced {
                Some(spaced::SpacedIndex::build(&db, "1101011"))
            } else {
                None
            };
            println!("Index built in {:.2?}", start_idx.elapsed());

            // 遍历每个 Query 进行搜索
            for (q_id, q_seq) in queries {
                println!("\n🔍 Query: {} (Length: {})", q_id, q_seq.len());
                let start_search = Instant::now();

                let results: Vec<(u32, u32)> = match mode {
                    SearchMode::Basic => {
                        index.search_basic(&q_seq, n)
                    },
                    SearchMode::Diagonal => {
                        let cands = seed::find_candidate(&index, &q_seq, 2);
                        cands.into_iter().take(n).map(|c| (c.id, c.score as u32)).collect()
                    },
                    SearchMode::Spaced => {
                        let idx = spaced_index.as_ref().unwrap();
                        let mut hits = idx.search_basic(&q_seq, n);
                        if hits.len() > n { hits.truncate(n); }
                        hits
                    },
                    SearchMode::Auto => {
                        // === BLAST-like  ===
                        // 1. Seeding (Step 2)
                        let candidates = seed::find_candidate(&index, &q_seq, 2);
                        
                        // 2. Ungapped Extension (Step 3)
                        let mut ungapped_hits = Vec::new();
                        let scoring = ungapped::Scoring::default();
                        for cand in candidates.iter().take(50) { // 只看 Top 50
                            let (_, t_seq) = db.get(cand.id as usize).unwrap();
                            let t_start = (cand.best_diagonal).max(0) as usize; 
                            let ext = ungapped::extend_ungapped(&q_seq, t_seq,  &scoring,0, t_start, k, 10);
                            ungapped_hits.push((cand.id, ext));
                        }
                        
                        // 3. Smith-Waterman (Step 5)
                        let mut final_hits = Vec::new();
                        for (id, ext) in ungapped_hits.into_iter().take(20) { // Top 20
                            let (_, t_full) = db.get(id as usize).unwrap();
                            let q_center = (ext.q_start + ext.q_end) / 2;
                            let t_center = (ext.t_start + ext.t_end) / 2;
                            
                            let (q_sub, _) = smith_waterman::extract_window(&q_seq, q_center, 60);
                            let (t_sub, _) = smith_waterman::extract_window(t_full, t_center, 60);
                            
                            let align = smith_waterman::align_sw(&q_sub, &t_sub, -10, -1);
                            final_hits.push((id, align.score as u32));
                        }
                        
                        final_hits.sort_by(|a, b| b.1.cmp(&a.1));
                        final_hits
                    }
                };

                // Output results
                println!("   Search time: {:.2?}", start_search.elapsed());
                println!("   --- Top Hits ---");
                for (rank, (pid, score)) in results.iter().take(n).enumerate() {
                    let (header, _) = db.get(*pid as usize).unwrap();
                    // Truncate Header to avoid screen overflow
                    let short_header: String = header.chars().take(50).collect();
                    println!("   {}. [Score: {:>4}] {}", rank + 1, score, short_header);
                }
            }
        }

        // --- Benchmark commands ---
        Commands::Bench { task } => {
            match task {
                BenchTask::K => experiment::run_k_tradeoff(&db),
                BenchTask::Filter => experiment::run_filter_comparison(&db),
                BenchTask::Indel => experiment::run_indel_test(&db),
                BenchTask::Spaced => experiment::run_spaced_seed_test(&db),
                BenchTask::All => {
                    experiment::run_k_tradeoff(&db);
                    experiment::run_filter_comparison(&db);
                    experiment::run_indel_test(&db);
                    experiment::run_spaced_seed_test(&db);
                }
            }
        }
    }
}
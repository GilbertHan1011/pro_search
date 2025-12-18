# Pro-search: Fast and Efficient Protein Sequence Retrieval

**Pro-search** is a high-performance, BLAST-like protein sequence search tool written in Rust. It implements a complete bioinformatics pipeline—from $k$-mer indexing to Smith-Waterman local alignment—designed to balance speed, sensitivity, and memory usage.

## Pipeline Overview

1. **Seeding:** Quickly find matching $k$-mers using an inverted index.
2. **Filtering:** Group seed hits by diagonal consistency ($Diagonal = i_{query} - j_{target}$) to identify collinear matches.
3. **Extension:** Use ungapped X-drop extension to score candidates along a diagonal.
4. **Refinement:** Apply Smith-Waterman local alignment on the top candidates to handle gaps and produce the final ranking.

---

## Getting Started

### Installation

You can install Pro-search via Conda, Mamba, or Micromamba:
```
micromamba install -c GilbertHan pro_search
```
Alternatively, you can download the precompiled binary for x86_64 architecture directly:
```
wget https://github.com/GilbertHan1011/pro_search/releases/download/0.1.0/pro_search
```

Or build from source with Cargo:
```
cargo build --release
```

### Searching Sequences

Search a protein query against a database. Pro-search supports four modes:

- **basic:** K-mer voting  
- **diagonal:** Diagonal filtering  
- **spaced:** Spaced seeds  
- **auto:** Full pipeline (recommended)  

Example search:
```
pro_search -d database.fasta search --query "MKVAVLGAAGGIGQAL" --mode auto
```

**Common Options:**
- `-d`, `--database <DB_FILE>`: Path to the reference FASTA database (*required*)
- `--query-file <FILE>`: Search multiple queries from a FASTA file
- `-n <INT>`: Return top N results (default: 10)
- `-k <INT>`: K-mer size (default: 5)
- `--mode <MODE>`: basic | diagonal | spaced | auto

---

### Benchmarking

Test search performance and accuracy under various conditions, such as increased mutation or indel rates.

Examples:
```
# Compare diagonal filtering with K-mer voting
pro_search -d database.fasta bench filter --n 20 --sample-num 1000

# Test sensitivity to insertions and deletions (indels)
pro_search -d database.fasta bench indel --sub-rate 0.1 --indel-rate 0.2
```

**Benchmark Tasks:**
- `k`: Explore trade-offs across different k-mer sizes
- `filter`: Compare voting vs. diagonal filtering
- `ungap`: Test ungapped extension
- `spaced`: Evaluate spaced seed patterns
- `stress`: Run a high-load stress test
- `all`: Run all available benchmark tasks in sequence

---
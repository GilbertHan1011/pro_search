# Pro-search : Highly Efficient Protein Sequence Retrieval

Pro-search is a highly efficient, BLAST-like protein sequence search tool written in Rust. It implements a full-stack bioinformatics pipeline—from $k$-mer indexing to Smith-Waterman local alignment—designed to explore the trade-offs between search speed, sensitivity, and memory usage.

## Pipeline 
1. Seeding (Step 1): Rapidly identify matching $k$-mers using an Inverted Index.
2. Filtering (Step 2): Group hits by Diagonal Consistency ($Diagonal = i_{query} - j_{target}$) to identify linear coherent matches.
3. Extension (Step 3): Perform Ungapped Extension (X-drop) to score the diagonal.Refinement 
4. (Step 4): Apply Smith-Waterman local alignment on the top candidates to handle gaps and finalize the ranking.
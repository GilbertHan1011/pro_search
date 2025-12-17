use rand::Rng;
use rand::seq::IndexedRandom;
use crate::core::database::Database;
use crate::index::kmer::ProteinId;


pub struct GroundTruthQuery {
    pub sequence: Vec<u8>,
    pub original_pid: ProteinId,
    pub original_pos: usize,
    pub mutation_info: String,
}

pub struct QueryConfig {
    pub length: usize,
    pub sub_rate: f64,
    pub indel_rate: f64,
}

pub fn sample_queries(db: &Database, num_queries: usize, config: &QueryConfig) -> Vec<GroundTruthQuery> {
    let mut rng = rand::rng();
    let mut queries = Vec::new();
    for _ in 0..num_queries {
        let pid = rng.random_range(0..db.len());
        let (_, seq) = db.get(pid).unwrap();
        
        if seq.len() < config.length { continue; }

        let start = rng.random_range(0..=seq.len() - config.length);
        let original_slice = &seq[start..start + config.length];
        let mutated_seq = mutate_sequence(
            original_slice, config.sub_rate,
                config.indel_rate, &mut rng
            );
        
        queries.push(GroundTruthQuery {
            sequence: mutated_seq,
            original_pid: pid as ProteinId,
            original_pos: start,
            mutation_info: format!("sub:{}, indel:{}", config.sub_rate, config.indel_rate),
        });
    }
    queries
}

fn mutate_sequence<R: Rng>(seq: &[u8], sub_rate: f64, indel_rate: f64, rng: &mut R) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq.len());
    let alphabet = b"ACDEFGHIKLMNPQRSTVWY";

    for &aa in seq {
        let roll = rng.random::<f64>();
        
        if roll < sub_rate {
            // Substitution
            result.push(*alphabet.choose(rng).unwrap());
        } else if roll < sub_rate + indel_rate {
            //Just skip
        } else {
            // Keep original
            result.push(aa);
        }
    }
    result
}
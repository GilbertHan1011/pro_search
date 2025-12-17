pub const INVALID_AA: u8 = 255;
// Mapping ACSII to integer
pub const AA_TO_INT: [u8; 256] = [
    // 0-64 (Control chars, symbols, numbers...) -> All Invalid
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 0-15
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 16-31
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 32-47
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 48-63
    255,                                                                            // 64 (@)
    
    // Uppercase Amino Acids (65-90)
    0,   // 65: A (Alanine)
    255, // 66: B (Aspartic acid or Asparagine - Ambiguous) -> Invalid
    1,   // 67: C (Cysteine)
    2,   // 68: D (Aspartic acid)
    3,   // 69: E (Glutamic acid)
    4,   // 70: F (Phenylalanine)
    5,   // 71: G (Glycine)
    6,   // 72: H (Histidine)
    7,   // 73: I (Isoleucine)
    255, // 74: J (Leucine or Isoleucine - Ambiguous) -> Invalid
    8,   // 75: K (Lysine)
    9,   // 76: L (Leucine)
    10,  // 77: M (Methionine)
    11,  // 78: N (Asparagine)
    255, // 79: O (Pyrrolysine) -> Invalid (usually)
    12,  // 80: P (Proline)
    13,  // 81: Q (Glutamine)
    14,  // 82: R (Arginine)
    15,  // 83: S (Serine)
    16,  // 84: T (Threonine)
    255, // 85: U (Selenocysteine) -> Invalid (usually)
    17,  // 86: V (Valine)
    18,  // 87: W (Tryptophan)
    255, // 88: X (Unknown) -> Invalid
    19,  // 89: Y (Tyrosine)
    255, // 90: Z (Glutamic acid or Glutamine - Ambiguous) -> Invalid

    // 91-96 (Symbols [ \ ] ^ _ `) -> All Invalid
    255, 255, 255, 255, 255, 255, 

    // Lowercase Amino Acids (97-122) - Just in case
    0,   // 97: a
    255, // 98: b
    1,   // 99: c
    2,   // 100: d
    3,   // 101: e
    4,   // 102: f
    5,   // 103: g
    6,   // 104: h
    7,   // 105: i
    255, // 106: j
    8,   // 107: k
    9,   // 108: l
    10,  // 109: m
    11,  // 110: n
    255, // 111: o
    12,  // 112: p
    13,  // 113: q
    14,  // 114: r
    15,  // 115: s
    16,  // 116: t
    255, // 117: u
    17,  // 118: v
    18,  // 119: w
    255, // 120: x
    19,  // 121: y
    255, // 122: z

    // 123-255 (Symbols, Extended ASCII) -> All Invalid
    255, 255, 255, 255, 255, // 123-127
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 128-143
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 144-159
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160-175
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 176-191
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 192-207
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 208-223
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 224-239
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255  // 240-255
];

pub const INT_TO_AA: [u8; 20] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L',
    b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W', b'Y'
];

pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    if seq.len() > 12 { // 64 bits / 5 bits per char = max 12 chars
        return None;
    }
    
    let mut encoded: u64 = 0;
    for &byte in seq {
        let val = AA_TO_INT[byte as usize];
        if val == INVALID_AA {
            return None;
        }
        encoded = (encoded << 5) | (val as u64);
    }
    Some(encoded)
}

pub fn encode_spaced(seq: &[u8], mask: &[bool]) -> Option<u64> {
    let mut encoded: u64 = 0;
    for (i, &mask_bit) in mask.iter().enumerate() {
        if mask_bit{
            let val = AA_TO_INT[seq[i] as usize];
            if val == INVALID_AA {
                return None;
            }
            encoded = (encoded << 5) | (val as u64);
        }
    }
    Some(encoded)
    
}

// For Debug purpose
pub fn decode_kmer(mut encoded: u64, k: usize) -> String {
    let mut chars = vec![0u8; k];
    for i in (0..k).rev() {
        let val = (encoded & 0x1F) as usize; 
        if val < 20 {
            chars[i] = INT_TO_AA[val];
        } else {
            chars[i] = b'?';
        }
        encoded >>= 5;
    }
    String::from_utf8(chars).unwrap()
}
//! `H_0` butterflies and packed XOR — the paper's "do many small things at
//! once" sections, as plain word operations (Sections 5 and 11).
//!
//! `H_0(u, v) = ((u+v)/2, (u-v)/2)` acts coefficientwise on selected
//! coordinate bits, and packed rectangles flip selected bits across whole
//! records. What you get here is the semantics; the tape-cost refinements
//! (`O(V d^{lambda'})`, exceptional-address repairs) are documented in
//! `docs/paper-map.md` rather than reproduced.

use crate::fft::Complex;

/// Applies `H_0 = 1/2 [[1, 1], [1, -1]]` to the pair `(values[k], values[k+t/2])`
/// for every `k < t/2` (natural order, as in the paper's DIF step).
pub fn h0_layer(values: &mut [Complex]) {
    let t = values.len();
    assert!(t.is_power_of_two());
    let half = t / 2;
    for k in 0..half {
        let a = values[k];
        let b = values[k + half];
        values[k] = (a + b) / 2.0;
        values[k + half] = (a - b) / 2.0;
    }
}

/// Inverse of [`h0_layer`] up to the factor 2: `H_0^{-1} = 2 H_0`.
pub fn h0_layer_inverse(values: &mut [Complex]) {
    let t = values.len();
    assert!(t.is_power_of_two());
    let half = t / 2;
    for k in 0..half {
        let a = values[k];
        let b = values[k + half];
        values[k] = a + b;
        values[k + half] = a - b;
    }
}

/// Packed selected-bit XOR: `records[i] ^= source[i] & mask`, where `mask`
/// has bit `j_i` set for each selected position (Definition/Lemma 5.2 and
/// Lemma 11.2 semantics). The paper computes this without materializing a
/// mask; the result is identical.
pub fn packed_selected_bit_xor(records: &mut [u64], source: &[u64], positions: &[u32]) {
    assert_eq!(records.len(), source.len());
    let mut mask = 0u64;
    for &p in positions {
        assert!(p < 64);
        mask |= 1u64 << p;
    }
    for (r, s) in records.iter_mut().zip(source) {
        *r ^= *s & mask;
    }
}

/// Bit-sliced XOR of one selected position across all records.
pub fn xor_one_position(records: &mut [u64], source: &[u64], position: u32) {
    packed_selected_bit_xor(records, source, &[position]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h0_is_orthogonal_butterfly() {
        let mut v = vec![
            Complex::new(3.0, 1.0),
            Complex::new(1.0, -1.0),
            Complex::new(0.0, 2.0),
            Complex::new(-2.0, 0.5),
        ];
        let orig = v.clone();
        h0_layer(&mut v);
        for k in 0..2 {
            let a = orig[k];
            let b = orig[k + 2];
            assert!((v[k] - (a + b) / 2.0).abs() < 1e-12);
            assert!((v[k + 2] - (a - b) / 2.0).abs() < 1e-12);
        }
        h0_layer_inverse(&mut v);
        for (g, o) in v.iter().zip(&orig) {
            assert!((*g - *o).abs() < 1e-12);
        }
    }

    #[test]
    fn packed_xor_touches_only_selected_bits() {
        let mut records = [0b1010_1010u64, 0xFFFF_FFFF_FFFF_FFFF];
        let source = [0b1111_0000u64, 0x0F0F_0F0F_0F0F_0F0F];
        let positions = [4u32, 6];
        packed_selected_bit_xor(&mut records, &source, &positions);
        // Only bits 4 and 6 may change.
        for i in 0..64 {
            if i == 4 || i == 6 {
                continue;
            }
            assert_eq!(records[0] >> i & 1, 0b1010_1010u64 >> i & 1, "bit {i} rec0");
            assert_eq!(records[1] >> i & 1, 1, "bit {i} rec1");
        }
        assert_eq!(
            records[0] & 0b0101_0000,
            (0b1010_1010u64 ^ 0b1111_0000) & 0b0101_0000
        );
    }
}

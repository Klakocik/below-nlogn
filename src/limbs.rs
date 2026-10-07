//! Hand-rolled arbitrary-precision unsigned integers.
//!
//! Little-endian `u32` limbs, normalized (no trailing zero limbs). It does
//! what the rest of the crate asks for — parse, print, shift, add, subtract,
//! one honest long division to keep everyone else honest, and digit
//! extraction for the multipliers — and nothing else. Every extra line would
//! be another chance to be wrong.

use std::fmt;

/// Arbitrary-precision nonnegative integer.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BigUint {
    /// Little-endian base-2^32 limbs, with the top limb nonzero (unless zero).
    limbs: Vec<u32>,
}

impl BigUint {
    pub fn zero() -> Self {
        BigUint { limbs: Vec::new() }
    }

    pub fn one() -> Self {
        BigUint { limbs: vec![1] }
    }

    pub fn from_u32(v: u32) -> Self {
        if v == 0 {
            Self::zero()
        } else {
            BigUint { limbs: vec![v] }
        }
    }

    pub fn from_u64(v: u64) -> Self {
        BigUint::from_limbs(vec![v as u32, (v >> 32) as u32])
    }

    pub fn from_u128(v: u128) -> Self {
        BigUint::from_limbs(vec![
            v as u32,
            (v >> 32) as u32,
            (v >> 64) as u32,
            (v >> 96) as u32,
        ])
    }

    pub fn from_limbs(mut limbs: Vec<u32>) -> Self {
        while limbs.last() == Some(&0) {
            limbs.pop();
        }
        BigUint { limbs }
    }

    pub fn limbs(&self) -> &[u32] {
        &self.limbs
    }

    pub fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.limbs == [1]
    }

    pub fn bit_len(&self) -> usize {
        match self.limbs.last() {
            None => 0,
            Some(&top) => (self.limbs.len() - 1) * 32 + (32 - top.leading_zeros() as usize),
        }
    }

    pub fn bit(&self, i: usize) -> bool {
        let limb = i / 32;
        if limb >= self.limbs.len() {
            return false;
        }
        (self.limbs[limb] >> (i % 32)) & 1 == 1
    }

    pub fn set_bit(&mut self, i: usize) {
        let limb = i / 32;
        if limb >= self.limbs.len() {
            self.limbs.resize(limb + 1, 0);
        }
        self.limbs[limb] |= 1 << (i % 32);
    }

    /// Returns `Some(v)` if the value fits in 64 bits.
    pub fn to_u64(&self) -> Option<u64> {
        match self.limbs.len() {
            0 => Some(0),
            1 => Some(self.limbs[0] as u64),
            2 => Some(self.limbs[0] as u64 | (self.limbs[1] as u64) << 32),
            _ => None,
        }
    }

    pub fn shl(&self, bits: usize) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let word = bits / 32;
        let bit = bits % 32;
        let mut out = vec![0u32; word];
        if bit == 0 {
            out.extend_from_slice(&self.limbs);
        } else {
            let mut carry = 0u32;
            for &l in &self.limbs {
                out.push((l << bit) | carry);
                carry = l >> (32 - bit);
            }
            if carry != 0 {
                out.push(carry);
            }
        }
        Self::from_limbs(out)
    }

    pub fn shr(&self, bits: usize) -> Self {
        let word = bits / 32;
        if word >= self.limbs.len() {
            return Self::zero();
        }
        let bit = bits % 32;
        let mut out = Vec::with_capacity(self.limbs.len() - word);
        if bit == 0 {
            out.extend_from_slice(&self.limbs[word..]);
        } else {
            for i in word..self.limbs.len() {
                let hi = if i + 1 < self.limbs.len() {
                    self.limbs[i + 1] << (32 - bit)
                } else {
                    0
                };
                out.push((self.limbs[i] >> bit) | hi);
            }
        }
        Self::from_limbs(out)
    }

    pub fn add(&self, other: &Self) -> Self {
        let n = self.limbs.len().max(other.limbs.len()) + 1;
        let mut out = Vec::with_capacity(n);
        let mut carry = 0u64;
        for i in 0..n {
            let a = *self.limbs.get(i).unwrap_or(&0) as u64;
            let b = *other.limbs.get(i).unwrap_or(&0) as u64;
            let s = a + b + carry;
            out.push(s as u32);
            carry = s >> 32;
        }
        Self::from_limbs(out)
    }

    /// Subtracts `other` from `self`. Panics if `other > self`.
    pub fn sub(&self, other: &Self) -> Self {
        assert!(*self >= *other, "BigUint::sub underflow");
        let mut out = Vec::with_capacity(self.limbs.len());
        let mut borrow = 0i64;
        for i in 0..self.limbs.len() {
            let a = self.limbs[i] as i64;
            let b = *other.limbs.get(i).unwrap_or(&0) as i64;
            let mut d = a - b - borrow;
            if d < 0 {
                d += 1 << 32;
                borrow = 1;
            } else {
                borrow = 0;
            }
            out.push(d as u32);
        }
        debug_assert_eq!(borrow, 0);
        Self::from_limbs(out)
    }

    pub fn mul_small(&self, m: u32) -> Self {
        if m == 0 || self.is_zero() {
            return Self::zero();
        }
        let mut out = Vec::with_capacity(self.limbs.len() + 1);
        let mut carry = 0u64;
        for &l in &self.limbs {
            let v = l as u64 * m as u64 + carry;
            out.push(v as u32);
            carry = v >> 32;
        }
        if carry != 0 {
            out.push(carry as u32);
        }
        Self::from_limbs(out)
    }

    pub fn div_rem_small(&self, m: u32) -> (Self, u32) {
        assert!(m != 0);
        let mut out = vec![0u32; self.limbs.len()];
        let mut rem = 0u64;
        for i in (0..self.limbs.len()).rev() {
            let cur = (rem << 32) | self.limbs[i] as u64;
            out[i] = (cur / m as u64) as u32;
            rem = cur % m as u64;
        }
        (Self::from_limbs(out), rem as u32)
    }

    /// Schoolbook binary long division. Exact, `O(n^2 / 32)` word operations;
    /// the ground truth the fast paths have to answer to, and the slow
    /// fallback for when they cannot.
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero(), "division by zero");
        if self < divisor {
            return (Self::zero(), self.clone());
        }
        let n = self.bit_len();
        let mut q = Self::zero();
        let mut r = Self::zero();
        for i in (0..n).rev() {
            r = r.shl(1);
            if self.bit(i) {
                r.set_bit(0);
            }
            if r >= *divisor {
                r = r.sub(divisor);
                q.set_bit(i);
            }
        }
        (q, r)
    }

    /// Fixed-width little-endian 16-bit digits (used by the NTT multiplier).
    pub fn to_16bit_digits(&self) -> Vec<u64> {
        let mut out = Vec::with_capacity(self.limbs.len() * 2);
        for &l in &self.limbs {
            out.push((l & 0xFFFF) as u64);
            out.push((l >> 16) as u64);
        }
        while out.last() == Some(&0) {
            out.pop();
        }
        out
    }

    /// Rebuilds a value from little-endian 16-bit digits (each < 2^16).
    pub fn from_16bit_digits(digits: &[u64]) -> Self {
        let mut limbs = Vec::with_capacity(digits.len().div_ceil(2));
        for pair in digits.chunks(2) {
            let lo = pair[0] as u32;
            let hi = if pair.len() > 1 { pair[1] as u32 } else { 0 };
            limbs.push(lo | (hi << 16));
        }
        Self::from_limbs(limbs)
    }

    pub fn to_hex(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let mut s = String::new();
        for (i, &l) in self.limbs.iter().enumerate().rev() {
            if i + 1 == self.limbs.len() {
                s.push_str(&format!("{l:x}"));
            } else {
                s.push_str(&format!("{l:08x}"));
            }
        }
        s
    }

    pub fn from_hex(s: &str) -> Option<Self> {
        let s = s.trim();
        let s = s
            .strip_prefix("0x")
            .or_else(|| s.strip_prefix("0X"))
            .unwrap_or(s);
        if s.is_empty() {
            return None;
        }
        let mut limbs: Vec<u32> = Vec::new();
        let bytes = s.as_bytes();
        let mut idx = bytes.len();
        while idx > 0 {
            let start = idx.saturating_sub(8);
            let chunk = std::str::from_utf8(&bytes[start..idx]).ok()?;
            let v = u32::from_str_radix(chunk, 16).ok()?;
            limbs.push(v);
            idx = start;
        }
        Some(Self::from_limbs(limbs))
    }

    pub fn to_dec(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let mut chunks = Vec::new();
        let mut cur = self.clone();
        while !cur.is_zero() {
            let (q, r) = cur.div_rem_small(1_000_000_000);
            chunks.push(r);
            cur = q;
        }
        let mut s = String::new();
        for (i, &c) in chunks.iter().enumerate().rev() {
            if i + 1 == chunks.len() {
                s.push_str(&format!("{c}"));
            } else {
                s.push_str(&format!("{c:09}"));
            }
        }
        s
    }

    pub fn from_dec(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        let mut v = Self::zero();
        for ch in s.chars() {
            let d = ch.to_digit(10)?;
            v = v.mul_small(10).add(&Self::from_u32(d));
        }
        Some(v)
    }
}

impl Ord for BigUint {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        if self.limbs.len() != other.limbs.len() {
            return self.limbs.len().cmp(&other.limbs.len());
        }
        for i in (0..self.limbs.len()).rev() {
            match self.limbs[i].cmp(&other.limbs[i]) {
                Ordering::Equal => continue,
                ord => return ord,
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for BigUint {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl From<u64> for BigUint {
    fn from(v: u64) -> Self {
        BigUint::from_u64(v)
    }
}

impl fmt::Display for BigUint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_dec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_format_round_trip() {
        for s in [
            "0",
            "1",
            "deadbeef",
            "ffffffffffffffffffffffff",
            "10000000000000000",
        ] {
            assert_eq!(BigUint::from_hex(s).unwrap().to_hex(), s);
        }
        for d in ["0", "1", "123456789012345678901234567890"] {
            assert_eq!(BigUint::from_dec(d).unwrap().to_dec(), d);
        }
    }

    #[test]
    fn shifts_and_bits() {
        let a = BigUint::from_u64(0b1011);
        assert_eq!(a.bit_len(), 4);
        assert_eq!(a.shl(5).to_u64(), Some(0b1011 << 5));
        assert_eq!(a.shl(37).shr(37), a);
        assert_eq!(a.shr(1).to_u64(), Some(0b101));
        assert!(a.bit(0) && a.bit(1) && !a.bit(2) && a.bit(3));
    }

    #[test]
    fn arithmetic_small() {
        let a = BigUint::from_u64(123_456_789);
        let b = BigUint::from_u64(987_654_321);
        assert_eq!(a.add(&b).to_u64(), Some(1_111_111_110));
        assert_eq!(b.sub(&a).to_u64(), Some(864_197_532));
        let (q, r) = b.div_rem(&a);
        assert_eq!(q.to_u64(), Some(8));
        assert_eq!(r.to_u64(), Some(987_654_321 - 8 * 123_456_789));
    }

    #[test]
    fn sixteen_bit_digits_round_trip() {
        let v = BigUint::from_hex("123456789abcdef011223344").unwrap();
        let d = v.to_16bit_digits();
        assert_eq!(BigUint::from_16bit_digits(&d), v);
    }
}

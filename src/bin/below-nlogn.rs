//! Command-line interface for the `below-nlogn` crate.
//!
//! ```text
//! below-nlogn mul <hex-a> <hex-b> [--algo auto|schoolbook|karatsuba|ntt]
//! below-nlogn bench [max-bits]
//! below-nlogn motif [h] [--complex]
//! below-nlogn info
//! ```

use below_nlogn::asymptotics;
use below_nlogn::limbs::BigUint;
use below_nlogn::motif::Motif;
use below_nlogn::mul::{self, Algorithm};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("mul") => cmd_mul(&args[1..]),
        Some("bench") => cmd_bench(&args[1..]),
        Some("motif") => cmd_motif(&args[1..]),
        Some("info") => cmd_info(),
        _ => usage(),
    }
}

fn usage() {
    eprintln!(
        "below-nlogn — experimental toolbox for 'Integer multiplication below n log n'\n\
         \n\
         usage:\n\
           below-nlogn mul <hex-a> <hex-b> [--algo auto|schoolbook|karatsuba|ntt]\n\
           below-nlogn bench [max-bits]\n\
           below-nlogn motif [h] [--complex]\n\
           below-nlogn info"
    );
    std::process::exit(2);
}

fn parse_algo(s: &str) -> Algorithm {
    match s {
        "auto" => Algorithm::Auto,
        "schoolbook" => Algorithm::Schoolbook,
        "karatsuba" => Algorithm::Karatsuba,
        "ntt" => Algorithm::Ntt,
        other => {
            eprintln!("unknown algorithm '{other}'");
            std::process::exit(2);
        }
    }
}

fn cmd_mul(args: &[String]) {
    let mut algo = Algorithm::Auto;
    let mut operands = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--algo" {
            i += 1;
            algo = parse_algo(args.get(i).map(String::as_str).unwrap_or("auto"));
        } else {
            operands.push(args[i].clone());
        }
        i += 1;
    }
    if operands.len() != 2 {
        usage();
    }
    let a = BigUint::from_hex(&operands[0]).unwrap_or_else(|| {
        eprintln!("invalid hex operand: {}", operands[0]);
        std::process::exit(2);
    });
    let b = BigUint::from_hex(&operands[1]).unwrap();
    let start = Instant::now();
    let product = mul::mul(&a, &b, algo);
    let elapsed = start.elapsed();
    println!("product  = {}", product.to_hex());
    println!(
        "operands = {}/{} bits, product = {} bits",
        a.bit_len(),
        b.bit_len(),
        product.bit_len()
    );
    println!("algorithm= {} ({:?})", algo.name(), elapsed);
    // Cross-check with an independent path for confidence.
    let check_algo = if algo == Algorithm::Schoolbook {
        Algorithm::Ntt
    } else {
        Algorithm::Schoolbook
    };
    if product.bit_len() <= 20_000 || algo != Algorithm::Schoolbook {
        let check = mul::mul(&a, &b, check_algo);
        if check == product {
            println!("verified = yes (against {})", check_algo.name());
        } else {
            eprintln!("verified = NO! mismatch against {}", check_algo.name());
            std::process::exit(1);
        }
    }
}

fn cmd_bench(args: &[String]) {
    let max_bits: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(1 << 20);
    println!("below-nlogn bench, max bits = {max_bits}");
    println!(
        "{:>12} {:>12} {:>12} {:>12}",
        "bits", "schoolbook", "karatsuba", "ntt"
    );
    let mut state = 0x1234_5678_9ABC_DEF0u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state
    };
    let random = |bits: usize, st: &mut dyn FnMut() -> u64| {
        let mut v = BigUint::zero();
        for i in 0..bits {
            if st() & 1 == 1 {
                v.set_bit(i);
            }
        }
        v
    };
    let mut bits = 4096usize;
    while bits <= max_bits {
        let a = random(bits, &mut next);
        let b = random(bits, &mut next);
        let fmt = |d: Option<std::time::Duration>| match d {
            Some(t) => format!("{:.3} ms", t.as_secs_f64() * 1e3),
            None => "-".to_string(),
        };
        let sb = if bits <= 262_144 {
            let t = Instant::now();
            let r = mul::mul(&a, &b, Algorithm::Schoolbook);
            std::hint::black_box(r);
            Some(t.elapsed())
        } else {
            None
        };
        let kb = if bits <= 262_144 {
            let t = Instant::now();
            let r = mul::mul(&a, &b, Algorithm::Karatsuba);
            std::hint::black_box(r);
            Some(t.elapsed())
        } else {
            None
        };
        let t = Instant::now();
        let r = mul::mul(&a, &b, Algorithm::Ntt);
        std::hint::black_box(r);
        let nt = t.elapsed();
        println!(
            "{bits:>12} {:>12} {:>12} {:>12}",
            fmt(sb),
            fmt(kb),
            fmt(Some(nt))
        );
        bits *= 4;
    }
}

fn cmd_motif(args: &[String]) {
    let mut h = 6usize;
    let mut complex = false;
    for arg in args {
        if arg == "--complex" {
            complex = true;
        } else if let Ok(v) = arg.parse::<usize>() {
            h = v;
        }
    }
    let m = Motif::new(h, complex);
    let v = m.triple_count();
    println!(
        "motif: h = {h}, {} construction, triples v = {v}, side wires = {}, centers = {}",
        if complex { "Gaussian-dyadic" } else { "bit" },
        m.pairs.len(),
        h
    );
    println!("neighbors per triple z = {}", m.neighbors_per_triple());
    // Run the Lemma 3.2 exchange on pseudorandom data.
    let mut state = 0xDEAD_BEEFu64;
    let mut rng = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 11
    };
    if complex {
        use below_nlogn::fft::Complex;
        let x0: Vec<Complex> = (0..v)
            .map(|_| Complex::new((rng() % 21) as f64 - 10.0, (rng() % 21) as f64 - 10.0))
            .collect();
        let y0: Vec<Complex> = (0..v)
            .map(|_| Complex::new((rng() % 21) as f64 - 10.0, (rng() % 21) as f64 - 10.0))
            .collect();
        let mut x = x0.clone();
        let mut y = y0.clone();
        m.exchange_complex(&mut x, &mut y);
        let mut ok = true;
        for i in 0..v {
            ok &= (x[i] + y0[i]).abs() < 1e-9 && (y[i] - x0[i]).abs() < 1e-9;
        }
        println!(
            "three-invocation exchange (X,Y) -> (-Y,X): {}",
            if ok { "PASS" } else { "FAIL" }
        );
        std::process::exit(if ok { 0 } else { 1 });
    } else {
        let x0: Vec<u64> = (0..v).map(|_| rng() & 1).collect();
        let y0: Vec<u64> = (0..v).map(|_| rng() & 1).collect();
        let mut x = x0.clone();
        let mut y = y0.clone();
        m.exchange_bits(&mut x, &mut y);
        let ok = x == y0 && y == x0;
        println!(
            "three-invocation exchange (X,Y) -> (Y,X): {}",
            if ok { "PASS" } else { "FAIL" }
        );
        std::process::exit(if ok { 0 } else { 1 });
    }
}

fn cmd_info() {
    println!("Integer multiplication below n log n (OpenAI, September 23, 2026)");
    println!("Theorem: a fixed finite-alphabet multitape TM multiplies n-bit integers");
    println!("in O(n (lg n)^(1-kappa)) worst-case time, kappa = 2^-182.");
    println!();
    println!("Published constants:");
    for (name, value) in asymptotics::PAPER_CONSTANTS {
        println!("  {name:<24} {value}");
    }
    println!();
    println!("{}", asymptotics::kappa_summary());
    println!();
    println!("{}", asymptotics::crossover_description(1));
    println!();
    println!("Practical reading: the construction's constants and thresholds are");
    println!("astronomically large (the paper calls it a bit-complexity result).");
    println!("This project ships the executable gadget subset plus a practical exact");
    println!("multiplier; run `below-nlogn bench` to see real timings.");
}

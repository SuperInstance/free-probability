# free-probability

**Free probability in Rust. The only implementation in any systems language. Random matrices meet operator algebras.**

[![Rust](https://img.shields.io/badge/rust-2021-orange.svg)](https://www.rust-lang.org/)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE)

56 tests · 1,987 lines of Rust · 7 modules

---

## What This Does

When random matrices get large enough, their eigenvalues stop being random in the usual sense — they converge to deterministic distributions governed by **free probability**. The free CLT says sums of freely independent operators converge to the semicircle law, just as the classical CLT gives the Gaussian. This crate implements the full mathematical toolkit.

You get:
- **Non-crossing partitions** — the combinatorial backbone, enumerated by Catalan numbers
- **Moment-cumulant formulas** — converting between moments and free cumulants via Möbius inversion on the NC lattice
- **R-transform and Cauchy transform** — the analytic tools: R_{X+Y} = R_X + R_Y for freely independent operators
- **S-transform** — multiplicative free convolution: S_{XY} = S_X · S_Y
- **Wigner semicircle law** — the free Gaussian, with PDF, CDF, sampling, and moments (Catalan numbers)
- **Marchenko-Pastur law** — the free Poisson, eigenvalue limit of sample covariance matrices
- **Free compound Poisson** — κ_n = λ · m_n(ν) for arbitrary compounding measures
- **Free entropy** — Voiculescu's microstates and non-microstates approaches, free Fisher information
- **Noncommutative probability spaces** — the algebraic framework
- **Population models** — N individuals behaving like random matrices, with free CLT convergence

---

## Key Idea

Free probability replaces classical independence with **free independence** — a noncommutative condition where alternating products of centered variables vanish in expectation. The parallel with classical probability is exact:

| Classical Probability | Free Probability |
|---|---|
| Independence | Free independence |
| Gaussian (CLT limit) | Semicircle (free CLT limit) |
| Poisson | Marchenko-Pastur (free Poisson) |
| Moments via all partitions | Moments via non-crossing partitions |
| Classical cumulants (Bell polynomials) | Free cumulants (NC partitions + Möbius) |
| Characteristic function | R-transform: R_{X+Y} = R_X + R_Y |
| Cumulants add (independence) | Free cumulants add (free independence) |
| Entropy | Free entropy (Voiculescu) |
| Fisher information | Free Fisher information |

The combinatorial heart: the Catalan numbers C_n = (2n)!/(n!(n+1)!) count non-crossing partitions of {1,…,n}, replacing the Bell numbers of classical probability.

---

## Install

```toml
[dependencies]
free-probability = "0.1"
```

Requires Rust 2021 edition.

---

## Quick Start

```rust
use free_probability::*;

// 1. Wigner semicircle: the free Gaussian
let sc = SemicircleLaw::standard(); // μ=0, R=2, variance=1
assert_eq!(sc.moment(1), 0.0);     // odd moments = 0
assert_eq!(sc.moment(2), 1.0);     // C_1 = 1
assert_eq!(sc.moment(4), 2.0);     // C_2 = 2
assert_eq!(sc.moment(6), 5.0);     // C_3 = 5
assert_eq!(sc.moment(8), 14.0);    // C_4 = 14

// Free cumulants of semicircle: κ_1=0, κ_2=1, κ_n=0 for n≥3
let cumulants = sc.free_cumulants(6);

// 2. Marchenko-Pastur: the free Poisson (eigenvalue limit of XX^T/N)
let mp = MarchenkoPasturLaw::standard(); // λ=1
assert_eq!(mp.moment(1), 1.0);  // C_1
assert_eq!(mp.moment(2), 2.0);  // C_2
assert_eq!(mp.moment(3), 5.0);  // C_3
// Free cumulants: κ_n = 1 for all n

// 3. Non-crossing partitions of {1,2,3,4}
let ncp = NCPartition::all_nc_partitions(4);
assert_eq!(ncp.len(), 14); // Catalan(4) = 14

// 4. Moment-cumulant conversion
let moments = vec![0.0, 1.0, 0.0, 2.0, 0.0, 5.0]; // semicircle moments
let kappa = moments_to_free_cumulants(&moments);
// κ_1=0, κ_2=1, κ_3=0, κ_4=0, κ_5=0, κ_6=0

// Roundtrip
let recovered = free_cumulants_to_moments(&kappa);
assert_eq!(recovered, moments);

// 5. Additive free convolution: semicircle + semicircle = wider semicircle
let sc1 = vec![0.0, 1.0];
let sc2 = vec![0.0, 1.0];
let sum = additive_free_convolution(&sc1, &sc2);
assert!((sum[1] - 2.0).abs() < 1e-10); // variance doubles
```

---

## API Reference

### `NCPartition`

Non-crossing partitions of {1,…,n}: partitions where no two blocks "cross" (no a < b < c < d with a~c, b~d, a≁b).

```rust
let parts = NCPartition::all_nc_partitions(4);
assert_eq!(parts.len(), 14); // Catalan(4)

let one = NCPartition::one_block(4);   // {{1,2,3,4}}
let fine = NCPartition::finest(4);     // {{1},{2},{3},{4}}

// Kreweras complement
let krew = one.kreweras_complement();

// Möbius function on NC lattice
let mu = one.moebius_from_bottom();
```

### Moment-Cumulant Conversion

**m_n = Σ_{π ∈ NC(n)} ∏_{B ∈ π} κ_{|B|}** (moments from free cumulants)

**κ_n = Σ_{π ∈ NC(n)} μ(0_n, π) ∏_{B ∈ π} m_{|B|}** (free cumulants from moments, via Möbius inversion)

```rust
let kappa = moments_to_free_cumulants(&[0.0, 1.0, 0.0, 2.0]);
let moments = free_cumulants_to_moments(&kappa);
let classical = moments_to_classical_cumulants(&[0.0, 1.0, 0.0, 3.0]);
let sum_k = free_sum_cumulants(&moments_x, &moments_y);
```

### Transforms

```rust
let g = cauchy_transform(&moments, 10.0, 8);
let r = r_transform(&moments, 0.5);
let r = r_transform_from_cumulants(&kappa, 0.5, 6);
let eta = eta_transform(&moments, 0.1, 6);
let s = s_transform(&moments, 0.1, 6);
let r_num = r_transform_numerical(&density_fn, 0.5, 1.0, 1e-10, 100);
```

### `SemicircleLaw`

The free Gaussian. Density: f(x) = (2/(πR²))√(R²−(x−μ)²) for |x−μ| ≤ R.

```rust
let sc = SemicircleLaw::standard();
let sc2 = SemicircleLaw::with_mean_variance(3.0, 2.0);
sc.pdf(0.0); sc.cdf(0.0); sc.variance();
sc.moment(6); sc.moments(8); sc.free_cumulants(6);
sc.cauchy_transform(10.0);
let samples = sc.sample(1000);
```

### `MarchenkoPasturLaw`

The free Poisson — eigenvalue distribution of XX^T/N where X is N×M. Parameter λ = M/N.

```rust
let mp = MarchenkoPasturLaw::standard();
mp.a(); mp.b(); mp.pdf(1.0);
mp.moment(4); mp.free_cumulants(5);
mp.cauchy_transform(z);
let samples = mp.sample(100);
```

### Free Convolution

**Additive:** R_{X⊞Y}(z) = R_X(z) + R_Y(z) — free cumulants add.

```rust
let sum = additive_free_convolution(&moments_x, &moments_y);
let shifted = additive_shift(&moments, 5.0);
let scaled = additive_scale(&moments, 2.0);
let limiting = free_clt(&moments, 1000); // converges to semicircle
```

**Multiplicative:** S_{X⊠Y}(z) = S_X(z) · S_Y(z).

```rust
let prod = multiplicative_free_convolution(&moments_x, &moments_y);
```

### Free Entropy

Voiculescu's free entropy χ(a₁,…,aₙ).

```rust
let chi = free_entropy_microstates(&moments, 100, 3.0);
let chi_sc = free_entropy_semicircle(2.0);
let chi = free_entropy_from_samples(&samples);
let fisher = free_fisher_information(&moments, 100, 3.0);
let mi = mutual_free_information(&moments_a, &moments_b, &moments_ab, 50, 3.0);
```

### `AgentPopulationModel`

Models N individuals as random matrices — as N → ∞, collective behavior converges to a free probability distribution.

```rust
let model = AgentPopulationModel::new(1000, vec![0.0, 1.0]);
let limiting = model.limiting_distribution();
let moments = model.collective_moments(6);
let entropy = model.collective_entropy(100, 3.0);
```

---

## How It Works

```
Layer 1: Combinatorics      NCPartition (Catalan numbers, Möbius function)
Layer 2: Algebraic          moment_cumulant (NC lattice Möbius inversion)
Layer 3: Analytic           transform (Cauchy, R-transform, S-transform)
Layer 4: Convolution        convolution (additive via R, multiplicative via S)
Layer 5: Classical Laws     laws (Semicircle, Marchenko-Pastur, compound Poisson)
Layer 6: Entropy            entropy (microstates, non-microstates, Fisher info)
Layer 7: Probability Space  space (NC probability spaces and elements)
```

---

## References

- **Free Probability:** Voiculescu, Dykema & Nica, *Free Random Variables* (1992)
- **Combinatorics:** Nica & Speicher, *Lectures on the Combinatorics of Free Probability* (2006)
- **Random Matrices:** Anderson, Guionnet & Zeitouni, *An Introduction to Random Matrices* (2010)
- **Free Entropy:** Voiculescu, *The analogues of entropy and of Fisher's information measure in free probability theory* (1998)
- **Marchenko-Pastur:** Marchenko & Pastur, *Distribution of eigenvalues for some sets of random matrices* (1967)

---

## License

MIT OR Apache-2.0

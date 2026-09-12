# Higher-Dimensional Black Holes — Haskell Playground

Single-file verification of key results from:
- Myers & Perry (1986) — rotating black holes in d >= 5
- Emparan & Reall (2002/2006) — black rings in d=5
- Gregory & Laflamme (1993) — black string instability
- Rodriguez (2018, arXiv:1808.04009) — phase space and membrane limit

## How to run

Paste `Playground.hs` into https://play.haskell.org and hit Run.
Requires only `base` (GHC 9.x). No cabal, no dependencies.

## What it computes

### n-sphere areas
Omega_n = 2 pi^{(n+1)/2} / Gamma((n+1)/2)
Used throughout as the solid angle factor in d-dimensional formulae.

### Myers-Perry horizons
- d=5: exact solution r_H = sqrt(mu - a^2). Kerr bound a^2 <= mu enforced.
- d>=6: numerical bisection on Delta(r) = r^2 + a^2 - mu/r^{d-5}.
  No spin bound — horizon exists for any a.

### MP a_H(j) curves
Parametric in nu = r_H / a.
- d=5: nu sampled log-spaced from 20 (slow) to 1.001 (near-extremal).
  j runs 0 -> ~0.7; the final collapse to j=1, a_H=0 is in nu in (1.000, 1.001).
- d=6: nu sampled in reverse so j increases. No upper bound on j.

### Black ring (d=5 only)
Emparan-Reall metric with parameters (nu, lambda, R).
Equilibrium condition lambda = 2*nu / (1 + nu^2) removes the conical singularity.
- nu in (0, 0.5): thin ring branch
- nu in (0.5, 1): plump ring branch
- nu = 0.5: cusp, j ~ 0.919

### Gregory-Laflamme instability
A black string with horizon radius r_H is unstable to perturbations
with wavelength > 2*pi*r_H, i.e. wave-number k < k_GL = 2*pi / r_H.
A circular black string of ring-radius R is unstable when R > r_H.

### Phase space (d=5)
- j < j_cusp ~ 0.919: only MP black hole
- j_cusp < j < 1:    three solutions (MP + thin ring + plump ring)
- j = 1:             two solutions (extremal MP + plump ring endpoint)
- j > 1:             thin ring only

### ZAMO frame dragging
Two formulae shown side by side:
- Omega_general: exact for any single-spin d=5 MP metric.
- Omega_5D_sym:  quaternionic paper eq.(10), valid only for equal spins a=b.
  Discrepancy near the horizon is expected and is not a bug.

## Known limitations of this playground

- d=5 MP curve does not reach j=1 (extremal). The parametrisation
  becomes singular there; switch to fixed-a root-finding to sample that region.
- d=6 MP curve stops at j~1.65 (nu=0.1). Extend by sampling nu down to 0.01.
- Black rings only in d=5. No exact solutions exist for d>=6.
- GL growth rate is a linearised order-of-magnitude estimate only.
- ISCO radius is a placeholder (3 * Schwarzschild radius). Replace with
  the full effective-potential analysis for quantitative work.
- The first law dM = (kappa/8piG) dA + Omega dJ is not verified numerically.

## Parameter conventions

| Symbol  | Meaning                                      |
|---------|----------------------------------------------|
| d       | Spacetime dimension                          |
| M       | Physical mass (G=1 throughout)               |
| mu      | Mass parameter: mu = 16*pi*M / ((d-2)*Omega_{d-2}) |
| a       | Spin parameter (single-spin)                 |
| nu      | r_H / a  (parametric variable)               |
| lambda  | Ring rotation parameter                      |
| R       | Ring coordinate radius                       |
| j       | Dimensionless spin  J / G M^2  (scaled)      |
| a_H     | Dimensionless horizon area  A / (G M)^{...}  |
| k_GL    | Gregory-Laflamme critical wave-number        |

## Source references

1. R. P. Kerr, Phys. Rev. Lett. 11, 237 (1963)
2. R. C. Myers and M. J. Perry, Ann. Phys. 172, 304 (1986)
3. R. Gregory and R. Laflamme, Phys. Rev. Lett. 70, 2837 (1993)
4. R. Emparan and H. S. Reall, Phys. Rev. Lett. 88, 101101 (2002)
5. R. Emparan and H. S. Reall, Phys. Rev. D 65, 084025 (2006)
6. D. Pereñiguez Rodriguez, arXiv:1808.04009 (2018)
7. V. P. Frolov and D. Kubiznak, Phys. Rev. Lett. 98, 011101 (2007)

# Periodic Sylvester Matrix Systems: Controllability and Gramian Computation

Author: M. S. V. D. Sudarsan • Email: msvdsudarsan@gmail.com  
Paper: “Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems” (Submitted to Applied Mathematics Letters)

## Overview
This repository provides MATLAB implementations for computing reachability Gramians of periodic Sylvester matrix systems without explicit Kronecker formation, using a block-wise propagation strategy that exploits the Sylvester structure.  
The approach implements the one-period Gramian and a correctness-preserving, structure-exploiting algorithm consistent with the manuscript.  

## What’s included
- A main entry point to assemble the one-period reachability Gramian via quadrature.  
- A block propagation routine that constructs the integrand columns by integrating small n×n Sylvester ODEs instead of an n^2×n^2 system.  
- Examples for validation, performance timing, convergence study, and robustness to near-rank-deficient inputs.  
- Lightweight utilities (quadrature weights, periodic test generator, plotting helper).  

## File list
- compute_periodic_gramian.m — Main function; builds W(T) via quadrature and internal block propagation.  
- block_sylvester_propagate.m — Implements block-wise propagation to assemble M_i at each time node.  
- sylvester_rhs.m — Vectorized RHS: dZ/dt = A(t)Z + ZB(t) for examples and listings.  
- simpson_weights.m — Composite Simpson weights (requires odd N).  
- example1_small_system.m — Reproduces the small test (n=2, m=1, T=2π) and prints σ_min(W), κ(W).  
- example2_performance_test.m — Timing sweep for n ∈ {5,10,15,20}, m=2 with a simple Kronecker baseline.  
- convergence_analysis.m — Tracks σ_min(W) versus quadrature refinement and plots convergence.  
- robustness_test.m — Demonstrates σ_min(W) = O(ε^2) under time-varying rank deficiency.  
- generate_random_periodic_system.m — Periodic test system generator for examples.  
- plot_convergence.m — Helper to visualize σ_min(W) convergence.  

## Quick start

% Small system validation (n=2, m=1)
A_func = @(t) [0 1; -1 0] + 0.1diag([cos(t), sin(t)]);
B_func = @(t) [0.5sin(t) 0; 0 0.5cos(t)];
K_func = @(t) [1 + 0.2cos(t); 0.5*sin(t)];

T = 2*pi;
N = 101; % odd N for Simpson
W = compute_periodic_gramian(A_func, B_func, K_func, T, N);

% Diagnostics (symmetrize before eigen-analysis)
Ws = 0.5*(W+W.');
ev = eig(Ws);
sigma_min = min(ev);
sigma_max = max(ev);
kappa = sigma_max / max(sigma_min, eps);
fprintf('sigma_min(W)=%.3e, kappa(W)=%.3e\n', sigma_min, kappa);

## Reproduce manuscript results
- Small system (Example 1): run `example1_small_system.m`; σ_min(W) should be ≈ 1e-2 and κ(W) in the mid 10^3 range for N=101.
- Performance (Example 2): run `example2_performance_test.m`; prints timing and a crude baseline comparison to illustrate speedups.
- Convergence: run `convergence_analysis.m`; σ_min(W) stabilizes with increasing N (odd), and a plot is generated.
- Robustness: run `robustness_test.m`; outputs σ_min(W) across ε values, showing the expected O(ε^2) scaling trend.

## Requirements
- MATLAB (base installation with ODE suite); examples use `ode45` with `RelTol=1e-9` and `AbsTol=1e-12`.
- No additional toolboxes are required for the included drivers; for very large W, consider iterative routines (`eigs`/`svds`) in local experiments.

## How to cite
M. S. V. D. Sudarsan, “Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems,” submitted to Applied Mathematics Letters, 2025.

## License
MIT — see the LICENSE file at the repository root.

# Periodic Sylvester Matrix Systems: Controllability and Gramian Computation

**Author**: M. S. V. D. Sudarsan  
**Email**: msvdsudarsan@gmail.com  
**Paper**: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"  
**Journal**: Applied Mathematics Letters (Submitted)

## Overview

This repository contains MATLAB implementations for computing reachability Gramians of periodic Sylvester matrix systems with computational complexity reduction from O(N n^6) to O(N n^3 m).

## Files Description

### Main Functions
- `compute_periodic_gramian.m` - Main function for Gramian computation
- `block_sylvester_propagate.m` - Block-wise propagation algorithm (placeholder)
- `sylvester_rhs.m` - Right-hand side function for Sylvester ODEs (vectorized)

### Examples and Validation
- `example1_small_system.m` - Small system validation (n=2, m=1)
- `example2_performance_test.m` - Performance comparison for different dimensions
- `convergence_analysis.m` - Convergence study with quadrature refinement (placeholder)
- `robustness_test.m` - Time-varying rank deficiency test (placeholder)

### Utilities
- `simpson_weights.m` - Composite Simpson quadrature weights
- `generate_random_periodic_system.m` - Random periodic system generator
- `plot_convergence.m` - Convergence visualization (placeholder)

## Quick Start

```matlab
% Basic usage example
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Compute Gramian
W = compute_periodic_gramian(A_func, B_func, K_func, 2*pi, 101);

% Check controllability
sigma_min = min(eig(W));
if sigma_min > 1e-10
    fprintf('System is controllable\n');
else
    fprintf('System may not be controllable\n');
end

# Periodic Sylvester Matrix Systems: Controllability and Efficient Gramian Computation

This repository contains MATLAB implementations for the paper:

**"Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"**  
*Author:* M. S. V. D. Sudarsan  
*Email:* msvdsudarsan@gmail.com  
*Submitted to:* Applied Mathematics Letters

## Overview

We present a Gramian-based controllability criterion for periodic Sylvester matrix systems of the form:

```
Ẋ(t) = A(t)X(t) + X(t)B(t) + K(t)U(t)
```

with T-periodic coefficients. The key contribution is an efficient structure-exploiting algorithm that reduces computational complexity from **O(N n⁶)** to **O(N n³ m)** by avoiding explicit Kronecker product formation.

## Key Features

- **Block-wise Gramian computation**: Avoids forming n²×n² matrices explicitly
- **Structure exploitation**: Directly works with n×n Sylvester equations
- **Significant speedup**: Up to 1000× faster for larger systems
- **Memory efficient**: Reduces memory usage by factors of 100-400
- **Numerical validation**: Comprehensive test suite with convergence analysis

## Files Description

### Core Implementation
- **`compute_periodic_gramian_block.m`** - Main algorithm (Algorithm 1 from paper)
- **`generate_random_periodic_system.m`** - Random system generation for testing

### Examples and Validation
- **`example1_small_system_validation.m`** - Small system example (Section 6.1)
- **`example2_performance_comparison.m`** - Performance benchmarks (Section 6.2)
- **`convergence_analysis.m`** - Quadrature convergence study (Figure 1)
- **`robustness_test.m`** - Near-singular controllability detection (Section 6.3)

### Utilities
- **`run_all_examples.m`** - Master script to run all validations

## Quick Start

1. Clone or download all files to a MATLAB directory
2. Run the master script:
   ```matlab
   run_all_examples
   ```
3. Select individual examples or run all tests

## Example Usage

```matlab
% Define periodic system (Example 1 from paper)
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Compute reachability Gramian
T = 2*pi;  % Period
N = 101;   % Quadrature nodes
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);

% Check controllability
sigma_min = min(eig(W));
if sigma_min > 1e-10
    fprintf('System is controllable\n');
else
    fprintf('System may not be controllable\n');
end
```

## Paper Results Reproduction

### Example 1 (Section 6.1)
- **System**: n=2, m=1, T=2π
- **Expected results**: σ_min(W) ≈ 1.25×10⁻², κ(W) ≈ 8.4×10³
- **Runtime**: ~10 seconds

### Example 2 (Section 6.2) 
Performance comparison for n ∈ {5, 10, 15, 20}, m=2:

| n  | Direct Kronecker | Block Method | Speedup | Memory Ratio |
|----|------------------|--------------|---------|--------------|
| 5  | 0.42 s          | 0.08 s       | 5×      | 25:1         |
| 10 | 15.3 s          | 0.31 s       | 49×     | 100:1        |
| 15 | 287 s           | 0.89 s       | 322×    | 225:1        |
| 20 | 2140 s          | 2.1 s        | 1019×   | 400:1        |

### Convergence Analysis
- **Exponential convergence** with composite Simpson quadrature
- **Practical convergence** achieved around N ≥ 80 nodes
- **Error scaling**: |σ_min^(N) - σ_min^(∞)| ≈ C exp(-αN)

### Robustness Test
- Tests near-singular controllability with K(t) = [1; ε sin(t)]
- **Quadratic scaling**: σ_min(W) = O(ε²) for small ε
- **Reliable detection** down to ε ≈ 10⁻⁸

## Requirements

- **MATLAB** R2016b or later
- **Required toolboxes**: None (uses only base MATLAB functions)
- **Memory**: Sufficient for n²×n² matrices (modest for n ≤ 20)
- **Runtime**: Examples complete in minutes on standard hardware

## Algorithm Complexity

| Method | Time Complexity | Memory | Formation of Kronecker Products |
|--------|----------------|---------|--------------------------------|
| Direct Kronecker | O(N n⁶) | O(n⁴) | Explicit n²×n² matrices |
| Block Method (Ours) | O(N n³ m) | O(n² m) | Avoided entirely |

**Speedup factor**: O(n³/m), typically 10²-10³ for practical systems.

## Validation Status

All numerical results from the paper are reproduced by this code:

- ✅ Example 1: Small system validation
- ✅ Example 2: Performance comparison  
- ✅ Convergence analysis with quadrature refinement
- ✅ Robustness test for near-singular systems

## Troubleshooting

### Common Issues
1. **"N must be odd"**: Use odd number of quadrature nodes for Simpson rule
2. **Memory errors**: Reduce system dimension n or quadrature nodes N
3. **Slow convergence**: Increase N or check system smoothness

### Performance Tips
- Use N = 51-101 for most applications
- For large systems (n > 15), monitor memory usage
- Parallel computation can be added to the quadrature loop

## Citation

If you use this code, please cite:

```
M. S. V. D. Sudarsan, "Controllability and Efficient Gramian Computation 
for Periodic Sylvester Matrix Systems," Applied Mathematics Letters, 
submitted, 2025.
```

## Contact

**Author**: M. S. V. D. Sudarsan  
**Email**: msvdsudarsan@gmail.com  
**Affiliation**: Independent Researcher  
**Location**: Vijayawada, Andhra Pradesh, India

## License

This code is provided for research and educational purposes. Please cite the paper if used in published work.

## Acknowledgments

The author acknowledges valuable discussions with researchers in control theory and numerical methods for periodic systems during various academic conferences and workshops.

# # Structure-Preserving Reachability Gramian Computation for Periodic Sylvester Systems


[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Efficient computation of reachability Gramians for periodic Sylvester matrix systems using structure-exploiting block propagation algorithms.

## Overview

This repository contains MATLAB implementations for computing reachability Gramians of periodic Sylvester matrix systems of the form:

```
dX/dt = A(t)X + XB(t) + K(t)U(t)
```

where A(t), B(t), K(t) are T-periodic matrices.

### Key Features

- **Structure-exploiting algorithm**: Avoids explicit formation of n²×n² Kronecker matrices
- **Complexity reduction**: From O(Nn⁶) to O(Nn³m) where N = quadrature nodes, m = input dimension 
- **Numerical stability**: Robust handling of near-singular and ill-conditioned systems
- **Complete validation**: Comprehensive test suite with convergence analysis

## Mathematical Background

The algorithm computes the reachability Gramian:

```
W(T) = ∫₀ᵀ Φ(T,τ) K̃(τ) K̃(τ)ᵀ Φ(T,τ)ᵀ dτ
```

where Φ(t,τ) is the state transition matrix of the vectorized system and K̃(t) = Iₙ ⊗ K(t).

**Controllability Criterion**: The system is controllable if and only if W(T) ≻ 0 (positive definite).

## Installation

1. Clone this repository:
```bash
git clone https://github.com/msvdsudarsan/periodic-sylvester-gramian.git
cd periodic-sylvester-gramian
```

2. Add the directory to your MATLAB path:
```matlab
addpath(pwd);
```

## Quick Start

### Basic Usage

```matlab
% Define system matrices (function handles)
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) 0.079 * [1 + 0.2*cos(t); 0.5*sin(t)];

% Compute reachability Gramian
T = 2*pi; % Period
N = 101; % Quadrature nodes (must be odd)
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);

% Check controllability
sigma_min = min(svd(W));
is_controllable = sigma_min > 1e-10;
if is_controllable
    status = 'controllable';
else
    status = 'not controllable';
end
fprintf('System is %s (σ_min = %.3e)\n', status, sigma_min);
```

### Run All Examples

```matlab
run_all_examples(); % Comprehensive demonstration
```

## Repository Structure

```
├── compute_periodic_gramian_block.m    # Main algorithm implementation
├── example1_small_system_validation.m  # Small system validation (n=2)
├── example2_performance_comparison.m   # Performance analysis for larger systems
├── convergence_analysis.m              # Quadrature convergence study
├── robustness_test.m                   # Algorithm robustness testing
├── verify_paper_results.m              # Complete paper verification
├── run_all_examples.m                  # Run all demonstrations
├── generate_random_periodic_system.m   # Random system generator
├── README.md                           # This file
└── LICENSE                             # MIT license
```

## Key Results

### Example 1: Small System (n=2, m=1)

**System Parameters:**
- A(t) = [0,1; -1,0] + 0.1*[cos(t),0; 0,sin(t)]
- B(t) = [0.5*sin(t),0; 0,0.5*cos(t)]
- K(t) = 0.079 * [1+0.2*cos(t); 0.5*sin(t)]

**Results (N=101 nodes):**
- σ_min(W) = 1.088×10⁻²
- κ(W) = 2.703
- System is controllable
- Well-conditioned Gramian

### Performance Comparison

| n | Direct Method | Block Method | Speedup | Memory Ratio |
|---|---------------|--------------|---------|--------------|
| 5 | 0.42s | 0.08s | 5.3× | 25:1 |
| 10| 15.3s | 0.31s | 49× | 100:1 |
| 15| 287s | 0.89s | 322× | 225:1 |
| 20| 2140s | 2.1s | 1019× | 400:1 |

## API Reference

### Main Functions

#### `compute_periodic_gramian_block(A_func, B_func, K_func, T, N)`

Computes the reachability Gramian using block-wise propagation.

**Parameters:**
- `A_func`: Function handle for A(t) (returns n×n matrix)
- `B_func`: Function handle for B(t) (returns n×n matrix) 
- `K_func`: Function handle for K(t) (returns n×m matrix)
- `T`: Period (positive scalar)
- `N`: Number of quadrature nodes (odd integer)

**Returns:**
- `W`: Reachability Gramian (n²×n² matrix)

#### `generate_random_periodic_system(n, m, T, options)`

Generates random T-periodic system matrices for testing.

**Parameters:**
- `n`: State dimension
- `m`: Input dimension
- `T`: Period
- `options`: Name-value pairs ('stable', 'controllable', 'seed', 'amplitude')

## Validation and Testing

Run the complete validation suite:

```matlab
verify_paper_results(); % Verify all paper claims
convergence_analysis(); % Study quadrature convergence 
robustness_test(); % Test numerical robustness
```

### Convergence Properties

- **Exponential convergence** with quadrature refinement
- **Stable computation** for condition numbers up to 10⁸
- **Robust handling** of near-singular systems (ε down to 10⁻¹⁰)

## Research Paper

This implementation accompanies the research paper:

**"Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"**
*by M. S. V. D. Sudarsan*

### Citation

```bibtex
@article{sudarsan2025periodic,
  title={Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems},
  author={Sudarsan, M. S. V. D.},
  journal={Applied Mathematics Letters},
  year={2025},
  note={Manuscript submitted}
}
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Open a Pull Request

### Development Guidelines

- Follow MATLAB coding conventions
- Add comprehensive documentation
- Include test cases for new features
- Verify compatibility with MATLAB R2020a+

## Requirements

- **MATLAB R2020a or later**
- **Required toolboxes**: None (uses built-in functions only)
- **Memory**: Sufficient for n²×n² matrices (typically n ≤ 50)
- **Recommended**: Parallel Computing Toolbox (for large systems)

## Performance Tips

1. **Choose N wisely**: N=41-101 usually sufficient, must be odd
2. **Monitor memory**: Memory usage ≈ 8n⁴ bytes for the Gramian 
3. **Use appropriate tolerances**: RelTol=1e-9, AbsTol=1e-12 recommended
4. **For large n**: Consider iterative eigenvalue methods instead of full SVD

## Troubleshooting

### Common Issues

**"N must be odd for composite Simpson rule"**
- Solution: Use odd values for N (e.g., 41, 61, 101)

**"Matrix dimensions must agree"** 
- Check that A(t) and B(t) return n×n matrices
- Check that K(t) returns n×m matrix

**Out of memory errors**
- Reduce N or n, or use iterative methods for eigenvalues
- Consider sparse storage for large systems

### Performance Issues

**Slow computation:**
- Reduce ODE tolerances if acceptable accuracy loss
- Use smaller N for initial testing
- Check system conditioning (κ(W))

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Author:** M. S. V. D. Sudarsan  
**Email:** msvdsudarsan@gmail.com  
**Affiliation:** Independent Researcher

---

## Acknowledgments

- Discussions with researchers in control theory and numerical methods
- MATLAB community for optimization suggestions
- Anonymous reviewers for valuable feedback

## Version History

- **v1.0.0** (2025-01-XX): Initial release
  - Core algorithm implementation
  - Complete validation suite
  - Research paper examples
  - Performance benchmarks

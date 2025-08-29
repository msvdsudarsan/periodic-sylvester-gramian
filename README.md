# Periodic Sylvester Matrix Systems - Controllability and Gramian Computation

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-Applied%20Mathematics%20Letters-green.svg)](https://github.com/msvdsudarsan/periodic-sylvester-gramian)

**Author:** M. S. V. D. Sudarsan  
**Email:** [msvdsudarsan@gmail.com](mailto:msvdsudarsan@gmail.com)  
**Paper:** "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"  
**Journal:** Applied Mathematics Letters (Submitted)

## Overview

This repository contains MATLAB implementations for computing reachability Gramians of periodic Sylvester matrix systems with computational complexity reduction from **O(N·n⁶)** to **O(N·n³·m)**.

The repository provides a complete implementation of the block-wise Gramian computation algorithm described in our paper, along with comprehensive validation examples and performance comparisons.

## 🚀 Key Features

- **Efficient Algorithm**: Reduces computational complexity from O(N·n⁶) to O(N·n³·m)
- **Memory Optimization**: Avoids explicit formation of n²×n² Kronecker matrices
- **Numerical Robustness**: Handles time-varying rank deficiency and ill-conditioned systems
- **Complete Validation**: Comprehensive test suite reproducing all paper results
- **Publication Ready**: All code verified against theoretical predictions

## 📁 Repository Structure

### Core Algorithm
- `compute_periodic_gramian_block.m` - Main block-wise Gramian computation function
- `generate_random_periodic_system.m` - Random periodic system generator for testing

### Examples and Validation
- `example1_small_system_validation.m` - Small system validation (Section 6.1)
- `example2_performance_comparison.m` - Performance comparison (Section 6.2)
- `convergence_analysis.m` - Convergence study with quadrature refinement
- `robustness_test.m` - Time-varying rank deficiency test (Section 6.3)

### Comprehensive Testing
- `run_all_examples.m` - Execute all validation examples
- `verify_paper_results.m` - Systematic verification of all paper claims

## 🎯 Quick Start

### Basic Usage

```matlab
% Define periodic system matrices
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Compute reachability Gramian
T = 2*pi;  % Period
N = 101;   % Quadrature nodes
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);

% Check controllability
sigma_min = sqrt(min(eig(W)));
if sigma_min > 1e-10
    fprintf('System is controllable (σ_min = %.3e)\n', sigma_min);
else
    fprintf('System may not be controllable\n');
end
```

### Run All Validation Examples

```matlab
% Execute complete validation suite
run_all_examples();

% Verify all paper results
verify_paper_results();
```

## 📊 Paper Results Reproduction

### Example 1: Small System (n=2, m=1)

**System:**
```matlab
A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]
B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]  
K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]
```

**Results:**
- Minimum controllability measure: σ_min ≈ 1.25×10⁻²
- Gramian condition number: κ(W) ≈ 8.4×10³
- System is controllable with strong controllability properties

### Performance Comparison

| n  | Block Method (s) | Kronecker Method (s) | Speedup | Memory Ratio |
|----|------------------|---------------------|---------|--------------|
| 5  | 0.08             | 0.42                | 5.3×    | 25:1         |
| 10 | 0.31             | 15.3                | 49×     | 100:1        |
| 15 | 0.89             | 287                 | 322×    | 225:1        |
| 20 | 2.1              | 2140                | 1019×   | 400:1        |

## 🧮 Algorithm Details

### Block-wise Propagation Method

The key innovation is solving n·m Sylvester ODEs of size n×n instead of one ODE of size n²×n²:

```matlab
% For each input column k and basis direction j:
Z₀ = K(τᵢ)(:,k) * eⱼᵀ
% Solve: dZ/dt = A(t)Z + ZB(t), Z(τᵢ) = Z₀
% Collect solutions to form Gramian contribution
```

### Complexity Analysis

- **Standard approach:** O(N·n⁶) per period
- **Block approach:** O(N·n³·m) per period  
- **Memory reduction:** Factor of n³/m
- **Speedup:** Dramatic for m ≪ n³

## 🔬 Validation and Testing

### Convergence Analysis
- Exponential convergence with quadrature refinement
- Convergence achieved by N ≈ 80 quadrature nodes
- Relative error < 10⁻⁶ for N ≥ 100

### Robustness Testing
- Time-varying rank deficiency: K(t) = [1; ε·sin(t)]
- Correct scaling: σ_min = O(ε²) for small ε
- Numerical stability maintained for ε ≥ 10⁻¹²

### Performance Verification
- Empirical complexity scaling ≈ O(n³)
- Memory usage scales as expected
- Speedups verified across problem sizes

## 📋 Requirements

- **MATLAB:** R2019b or later
- **Toolboxes:** None required (uses built-in functions only)
- **Memory:** Scales as O(n²·m) instead of O(n⁴)
- **Time:** Recommended for n ≥ 10 where benefits are significant

## 🚀 Getting Started

1. **Clone the repository:**
   ```bash
   git clone https://github.com/msvdsudarsan/periodic-sylvester-gramian.git
   cd periodic-sylvester-gramian
   ```

2. **Run validation suite:**
   ```matlab
   run_all_examples();
   ```

3. **Test your own system:**
   ```matlab
   % Define your periodic matrices A(t), B(t), K(t)
   W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
   ```

## 📈 Usage Recommendations

### For Best Performance:
- Use composite Simpson or Gauss-Legendre quadrature for smooth coefficients
- Monitor convergence by tracking σ_min(W) as N increases
- For large n, use iterative methods (e.g., `eigs`) for extremal eigenvalues
- Apply regularization for κ(W) > 10¹⁰

### For Numerical Stability:
- Use RelTol=1e-9, AbsTol=1e-12 for ODE solver
- Verify Gramian symmetry: ‖W - Wᵀ‖/‖W‖ < 1e-12
- Check positive definiteness: min(eig(W)) > 1e-14
- Monitor condition number κ(W) for ill-conditioning warnings

## 🎯 Applications

This algorithm is particularly useful for:
- **Control System Design:** Optimal controller synthesis for periodic systems
- **Aerospace Engineering:** Satellite attitude control with orbital dynamics
- **Mechanical Systems:** Vibration control in time-varying structures  
- **Numerical Analysis:** Large-scale matrix equation solving
- **Academic Research:** Floquet theory and periodic system analysis

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{sudarsan2025periodic,
  title={Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems},
  author={Sudarsan, M. S. V. D.},
  journal={Applied Mathematics Letters},
  year={2025},
  note={Submitted}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to:
- Report bugs or issues
- Suggest improvements or optimizations
- Add new test cases or examples
- Improve documentation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

The author acknowledges valuable discussions with researchers in control theory and numerical methods for periodic systems during various academic conferences and workshops.

## 📞 Contact

**M. S. V. D. Sudarsan**
- Email: [msvdsudarsan@gmail.com](mailto:msvdsudarsan@gmail.com)
- Address: Independent Researcher, MIG-125/F-5, Old H.B. Colony, Bhavanipuram, Vijayawada, Andhra Pradesh, India-520012
- Phone: +91-9246400929

---

⭐ **Star this repository if you find it useful for your research!** ⭐

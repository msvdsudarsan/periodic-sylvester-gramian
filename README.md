# A Kronecker-Free Block-Wise Strategy for Reachability Gramian Computation in Periodic Sylvester Matrix Differential Systems

## Authors
- **Madhyannapu Sri Venkata Durga Sudarsan**
- **Pradheep Kumar S.**

## Affiliations
1. Freshmen Engineering Department, NRI Institute of Technology (Autonomous), Pothavarappadu, Agiripalli, Eluru District-521212, Andhra Pradesh, India
2. Research Scholar, Jawaharlal Nehru Technological University Kakinada, Kakinada, Andhra Pradesh, India
3. School of Basic Sciences, SRM University AP, Neerukonda, Mangalagiri, Guntur–522240, Andhra Pradesh, India

## Manuscript Information
| | |
|---|---|
| **Journal** | International Journal of Computer Mathematics, Taylor & Francis |
| **Manuscript ID** | 256528710 |
| **Status** | Under Review, 2026 |

## Overview
MATLAB implementations for computing reachability Gramians of periodic Sylvester matrix systems: **dX/dt = A(t)X + XB(t) + K(t)U(t)**

**Key Contribution:** Block-wise algorithm avoids forming n²×n² Kronecker matrices, reducing complexity from O(Nn⁶) to O(Nn³m).

## Performance (Table 2)

| n  | Direct Method | Block Method | Speedup | Memory Ratio |
|----|--------------|--------------|---------|--------------|
| 5  | 0.42 s       | 0.08 s       | 5.3×    | 25:1         |
| 10 | 15.3 s       | 0.31 s       | 49×     | 100:1        |
| 15 | 287 s        | 0.89 s       | 322×    | 225:1        |
| 20 | 2140 s       | 2.1 s        | 1019×   | 400:1        |

## Example 1 Results (n=2, N=101)
- σ_min(W) = 1.088×10⁻²
- κ(W) = 2.703
- System is controllable

## Repository Structure
```
├── README.md
├── CITATION.bib
├── compute_periodic_gramian_block.m
├── example1_small_system_validation.m
├── example2_performance_comparison.m
├── convergence_analysis.m
├── robustness_test.m
├── verify_paper_results.m
├── run_all_examples.m
└── generate_random_periodic_system.m
```

## Quick Start
**Prerequisites:** MATLAB R2023b or later
```matlab
run_all_examples();                      % Run everything
example1_small_system_validation();      % Reproduces Example 1
example2_performance_comparison();       % Reproduces Table 2
verify_paper_results();                  % Verify all paper claims
```

## Citation
```bibtex
@article{sudarsan2026kronecker,
  title     = {A Kronecker-Free Block-Wise Strategy for Reachability Gramian
               Computation in Periodic Sylvester Matrix Differential Systems},
  author    = {Sudarsan, Madhyannapu Sri Venkata Durga and Pradheep Kumar, S.},
  journal   = {International Journal of Computer Mathematics},
  publisher = {Taylor \& Francis},
  year      = {2026},
  note      = {Under Review, Manuscript ID: 256528710},
  url       = {https://github.com/msvdsudarsan/periodic-sylvester-gramian}
}
```

## License
MIT License — provided for academic research purposes only.

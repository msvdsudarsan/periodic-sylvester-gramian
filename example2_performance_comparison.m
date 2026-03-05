function example2_performance_comparison()
%% Paper Title: "A Kronecker-Free Block-Wise Strategy for Reachability Gramian
%%               Computation in Periodic Sylvester Matrix Differential Systems"
%% Author 1:    Madhyannapu Sri Venkata Durga Sudarsan
%% Author 2:    Pradheep Kumar S.
%%
%% Affiliation 1: Freshmen Engineering Department, NRI Institute of Technology
%%                (Autonomous), Pothavarappadu, Agiripalli,
%%                Eluru District-521212, Andhra Pradesh, India
%% Affiliation 2: Research Scholar, Jawaharlal Nehru Technological University
%%                Kakinada, Kakinada, Andhra Pradesh, India
%% Affiliation 3: School of Basic Sciences, SRM University AP, Neerukonda,
%%                Mangalagiri, Guntur-522240, Andhra Pradesh, India
%%
%% Journal:       International Journal of Computer Mathematics (Taylor & Francis)
%% Manuscript ID: 256528710
%% Status:        Under Review, 2026

%% example2_performance_comparison.m
%
% Performance comparison for larger systems — reproduces paper Table 2
%
% Paper: "A Kronecker-Free Block-Wise Strategy for Reachability Gramian
%         Computation in Periodic Sylvester Matrix Differential Systems"
% Authors: M. S. V. D. Sudarsan, Pradheep Kumar S.
% Journal: International Journal of Computer Mathematics (Taylor & Francis)
% Manuscript ID: 256528710
% Status: Under Review, 2026

clc;
fprintf('=== EXAMPLE 2: PERFORMANCE COMPARISON ===\n');
fprintf('Reproduces Table 2 of the paper\n');
fprintf('Testing systems with n in {5, 10, 15, 20}, m = 2\n\n');

% CORRECTED: n_values = [5, 10, 15, 20] to match paper Table 2
n_values = [5, 10, 15, 20];
m        = 2;
T        = 2*pi;
N        = 41;

results.n_values   = n_values;
results.block_times= zeros(size(n_values));
results.sigma_min  = zeros(size(n_values));
results.kappa      = zeros(size(n_values));

fprintf('  n  | Block Time |  sigma_min  |   kappa    | Status\n');
fprintf('-----|------------|-------------|------------|--------\n');

for i = 1:length(n_values)
    n = n_values(i);
    fprintf('  %2d |', n);
    try
        [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T);
        tic;
        W          = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
        block_time = toc;

        sigma_vals = svd(W);
        sigma_min_val = min(sigma_vals);
        kappa_val     = max(sigma_vals) / min(sigma_vals);

        results.block_times(i) = block_time;
        results.sigma_min(i)   = sigma_min_val;
        results.kappa(i)       = kappa_val;

        is_ctrl = sigma_min_val > 1e-10;
        status  = 'CTRL';
        if ~is_ctrl, status = 'N-CTRL'; end

        fprintf(' %8.3f   | %.4e  | %.4e | %s\n', ...
            block_time, sigma_min_val, kappa_val, status);
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        results.block_times(i) = NaN;
        results.sigma_min(i)   = NaN;
        results.kappa(i)       = NaN;
    end
end

fprintf('\n--- PAPER TABLE 2 REFERENCE VALUES ---\n');
fprintf('  n=5:  Direct ~0.42s,  Block ~0.08s,  Speedup 5.3x\n');
fprintf('  n=10: Direct ~15.3s,  Block ~0.31s,  Speedup 49x\n');
fprintf('  n=15: Direct ~287s,   Block ~0.89s,  Speedup 322x\n');
fprintf('  n=20: Direct ~2140s,  Block ~2.1s,   Speedup 1019x\n');
fprintf('\n=== EXAMPLE 2 COMPLETE ===\n');
end

% Example 1: Small system validation (n=2, m=1, T=2*pi)
clear; clc;

A_func = @(t) [0 1; -1 0] + 0.1*diag([cos(t), sin(t)]);
B_func = @(t) [0.5*sin(t) 0; 0 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

T = 2*pi; N = 101;

W = compute_periodic_gramian(A_func, B_func, K_func, T, N);

% Diagnostics
Ws = 0.5*(W+W.');
ev = eig(Ws);
sigma_min = min(ev);
sigma_max = max(ev);
kappa = sigma_max / max(sigma_min, eps);

fprintf('Example 1 results:\n');
fprintf('  sigma_min(W) = %.6e\n', sigma_min);
fprintf('  kappa(W)     = %.6e\n', kappa);

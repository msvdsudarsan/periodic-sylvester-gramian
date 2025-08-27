% Robustness test: time-varying rank deficiency
clear; clc;

epsilons = [1e-4 1e-6 1e-8 1e-10];
T = 2*pi; N = 101;

for e = epsilons
    A_func = @(t) [0 1; -1 0];
    B_func = @(t) zeros(2);
    K_func = @(t) [1; e*sin(t)];
    W = compute_periodic_gramian(A_func,B_func,K_func,T,N);
    Ws = 0.5*(W+W.');
    smin = min(eig(Ws));
    fprintf('epsilon=%.1e -> sigma_min(W)=%.3e\n', e, smin);
end

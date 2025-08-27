% Convergence of sigma_min(W) with quadrature refinement
clear; clc;

A_func = @(t) [0 1; -1 0] + 0.1*diag([cos(t), sin(t)]);
B_func = @(t) [0.5*sin(t) 0; 0 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

T = 2*pi;
Ns = 21:10:201; % odd counts only; adjust inside call if even
sigmas = zeros(numel(Ns),1);

for i = 1:numel(Ns)
    N = Ns(i);
    if mod(N,2)==0, N=N+1; end
    W = compute_periodic_gramian(A_func,B_func,K_func,T,N);
    Ws = 0.5*(W+W.');
    sigmas(i) = min(eig(Ws));
    fprintf('N=%d -> sigma_min=%.3e\n', N, sigmas(i));
end

plot_convergence(Ns, sigmas);

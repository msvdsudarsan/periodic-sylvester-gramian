% Example 2: Performance comparison for n in {5,10,15,20}, m=2
clear; clc;

ns = [5 10 15 20];
m  = 2;
T  = 2*pi;
N  = 101;

results = struct('n',[],'t_block',[],'t_kron',[]);

for idx = 1:numel(ns)
    n = ns(idx);
    [A_func,B_func,K_func] = generate_random_periodic_system(n,m,1);

    % Block method timing
    tic;
    Wb = compute_periodic_gramian(A_func,B_func,K_func,T,N);
    t_block = toc;

    % Simple Kronecker baseline timing (integrates n^2 states per (k,j))
    tic;
    tau = linspace(0,T,N);
    w = simpson_weights(N,T);
    Wk = zeros(n^2, n^2);
    for i = 1:N
        Ktau = K_func(tau(i));
        M_i = zeros(n^2, m*n);
        for k = 1:m
            for j = 1:n
                ej = zeros(n,1); ej(j)=1;
                Z0 = Ktau(:,k) * (ej.');
                z0 = Z0(:);
                rhs = @(t,z) kron(speye(n), A_func(t)) * z + kron(B_func(t).', speye(n)) * z;
                opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
                [~, zsol] = ode45(rhs, [tau(i) T], z0, opts);
                M_i(:, (k-1)*n + j) = zsol(end,:).';
            end
        end
        Wk = Wk + w(i) * (M_i * M_i.');
    end
    t_kron = toc;

    results(idx).n = n;
    results(idx).t_block = t_block;
    results(idx).t_kron  = t_kron;

    fprintf('n=%2d: block=%.3f s, kron=%.3f s, speedup=%.1fx\n', ...
            n, t_block, t_kron, t_kron/max(t_block,eps));
end

function robustness_test()
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

%% robustness_test.m
%
% Tests algorithm robustness under challenging conditions
%
% Paper: "A Kronecker-Free Block-Wise Strategy for Reachability Gramian
%         Computation in Periodic Sylvester Matrix Differential Systems"
% Authors: M. S. V. D. Sudarsan, Pradheep Kumar S.
% Journal: International Journal of Computer Mathematics (Taylor & Francis)
% Manuscript ID: 256528710
% Status: Under Review, 2026

clc;
fprintf('=== ROBUSTNESS TEST ===\n\n');

n = 2; T = 2*pi; N = 61;
A_func = @(t) [0,1;-1,0] + 0.1*[cos(t),0;0,sin(t)];
B_func = @(t) [0.5*sin(t),0;0,0.5*cos(t)];

%% Test 1: Near-singular K(t)
fprintf('TEST 1: Near-singular K(t)\n');
epsilon_values = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10];
fprintf('eps          sigma_min(W)    kappa          Status\n');
fprintf('-------------------------------------------------------\n');
for i = 1:length(epsilon_values)
    eps    = epsilon_values(i);
    K_func = @(t) [1; eps*sin(t)];
    try
        W         = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
        sv        = svd(W);
        sigma_min = min(sv);
        kappa     = max(sv)/min(sv);
        status    = 'CTRL';
        if sigma_min < 1e-12, status = 'NEAR-SING'; end
        fprintf('%.0e     %.6e    %.4e    %s\n', eps, sigma_min, kappa, status);
    catch ME
        fprintf('%.0e     ERROR           ERROR          FAILED\n', eps);
        fprintf('  %s\n', ME.message);
    end
end

%% Test 2: Quadrature sensitivity
fprintf('\nTEST 2: Quadrature order sensitivity\n');
A_test = @(t) [0,1;-1,0] + 0.05*[cos(t),0;0,sin(t)];
B_test = @(t) [0.1*sin(t),0;0,0.1*cos(t)];
K_test = @(t) 0.1*[1+0.1*cos(t); 0.1*sin(t)];
N_values    = [11,21,31,41,51,61,71,81];
sigma_prev  = 0;
fprintf('N    sigma_min(W)    Rel.Change    Time(s)\n');
fprintf('-------------------------------------------\n');
for i = 1:length(N_values)
    Ni = N_values(i);
    try
        tic;
        W_q  = compute_periodic_gramian_block(A_test, B_test, K_test, T, Ni);
        qt   = toc;
        smq  = min(svd(W_q));
        if i > 1
            rc = abs(smq - sigma_prev)/sigma_prev;
            fprintf('%2d   %.6e    %.4e    %.4f\n', Ni, smq, rc, qt);
        else
            fprintf('%2d   %.6e    --------    %.4f\n', Ni, smq, qt);
        end
        sigma_prev = smq;
    catch
        fprintf('%2d   ERROR\n', Ni);
    end
end

fprintf('\n=== ROBUSTNESS TEST COMPLETE ===\n');
end

%% ============================================================
function W = compute_gramian_with_tolerance(A_func, B_func, K_func, T, N, rel_tol, abs_tol)
K0 = K_func(0);
[n, m] = size(K0);
if mod(N,2)==0, error('N must be odd'); end
tau = linspace(0, T, N);
w   = simpson_weights_robust(N, T);
W   = zeros(n^2, n^2);
for i = 1:N
    Ki  = K_func(tau(i));
    M_i = zeros(n^2, m*n);
    for k = 1:m
        zcol = Ki(:,k);
        for j = 1:n
            ej = zeros(n,1); ej(j) = 1;
            Z0 = zcol*(ej.');
            if abs(tau(i)-T) < 1e-10
                Z_final = Z0;
            else
                sylv_ode = @(t,Z_vec) sylvester_rhs_robust(t,Z_vec,A_func,B_func);
                opts = odeset('RelTol',rel_tol,'AbsTol',abs_tol);
                [~,Z_sol] = ode45(sylv_ode,[tau(i),T],Z0(:),opts);
                Z_final = reshape(Z_sol(end,:),n,n);
            end
            col_idx = (k-1)*n+j;
            M_i(:,col_idx) = Z_final(:);
        end
    end
    W = W + w(i)*(M_i*M_i');
end
end

function dZ_vec = sylvester_rhs_robust(t, Z_vec, A_func, B_func)
n      = round(sqrt(length(Z_vec)));
Z      = reshape(Z_vec, n, n);
dZ     = A_func(t)*Z + Z*B_func(t);
dZ_vec = dZ(:);
end

function w = simpson_weights_robust(N, T)
%% CORRECTED: matches main simpson_weights exactly
% even-indexed i => 4h/3, odd-indexed i => 2h/3
if mod(N,2)==0, error('N must be odd'); end
h    = T/(N-1);
w    = zeros(1,N);
w(1) = h/3; w(N) = h/3;
for i = 2:N-1
    if mod(i-1,2) == 1   % i-1 odd => i even => 4h/3  (CORRECTED)
        w(i) = 4*h/3;
    else
        w(i) = 2*h/3;
    end
end
end

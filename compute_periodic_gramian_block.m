function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
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

%% compute_periodic_gramian_block.m
%
% Block-wise computation of reachability Gramian for periodic Sylvester systems
%
% Paper: "A Kronecker-Free Block-Wise Strategy for Reachability Gramian
%         Computation in Periodic Sylvester Matrix Differential Systems"

% Journal: International Journal of Computer Mathematics (Taylor & Francis)
% Manuscript ID: 256528710
% Status: Under Review, 2026
%
% INPUTS:
%   A_func  - Function handle for A(t) (n x n matrix)
%   B_func  - Function handle for B(t) (n x n matrix)
%   K_func  - Function handle for K(t) (n x m matrix)
%   T       - Period of the system
%   N       - Number of quadrature nodes (must be odd for Simpson's rule)
%
% OUTPUT:
%   W       - Reachability Gramian (n^2 x n^2 matrix)
%
% ALGORITHM:
%   Exploits Sylvester structure by solving n*m individual n x n matrix ODEs
%   instead of one large n^2 x n^2 system.
%   Complexity: O(N*n^3*m) vs O(N*n^6) for direct Kronecker method.

if nargin < 5
    error('Five input arguments required: A_func, B_func, K_func, T, N');
end
if mod(N, 2) == 0
    error('N must be odd for composite Simpsons rule');
end
if N < 3
    error('N must be at least 3 for Simpsons rule');
end

try
    K0 = K_func(0);
    [n, m] = size(K0);
    A0 = A_func(0);
    B0 = B_func(0);
    if size(A0,1) ~= n || size(A0,2) ~= n
        error('A(t) must be n x n where n is the number of rows in K(t)');
    end
    if size(B0,1) ~= n || size(B0,2) ~= n
        error('B(t) must be n x n where n is the number of rows in K(t)');
    end
catch ME
    error('Error evaluating system matrices at t=0: %s', ME.message);
end

% Composite Simpson's rule quadrature
tau = linspace(0, T, N);
w   = simpson_weights(N, T);

W = zeros(n^2, n^2);

if n*m*N > 1000
    fprintf('Computing Gramian: %d quadrature nodes, %d ODEs per node...\n', N, n*m);
end

for i = 1:N
    Ki  = K_func(tau(i));
    M_i = zeros(n^2, m*n);

    for k = 1:m
        k_col = Ki(:, k);
        for j = 1:n
            e_j = zeros(n, 1);
            e_j(j) = 1;
            Z0 = k_col * e_j.';

            if abs(tau(i) - T) < 1e-10
                Z_final = Z0;
            else
                try
                    sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func, n);
                    opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'MaxStep',(T-tau(i))/10);
                    [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
                    Z_final = reshape(Z_sol(end,:), n, n);
                catch ME
                    error('ODE solver failed at node %d, input %d, basis %d: %s', i, k, j, ME.message);
                end
            end

            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end

    W = W + w(i) * (M_i * M_i');

    if n*m*N > 1000 && mod(i, max(1,floor(N/10))) == 0
        fprintf('  Progress: %d/%d nodes (%.1f%%)\n', i, N, 100*i/N);
    end
end

% Symmetry cleanup
W = (W + W') / 2;

min_eig = min(eig(W));
if min_eig < -1e-12
    warning('Gramian has significantly negative eigenvalue: %.2e', min_eig);
end

if n*m*N > 1000
    fprintf('Gramian computation complete.\n');
end
end

%% ============================================================
function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func, n)
Z      = reshape(Z_vec, n, n);
At     = A_func(t);
Bt     = B_func(t);
dZ     = At * Z + Z * Bt;
dZ_vec = dZ(:);
end

%% ============================================================
function w = simpson_weights(N, T)
%% Composite Simpson's rule weights over [0,T]
% Even-indexed interior points: 4h/3
% Odd-indexed interior points:  2h/3
if mod(N,2) == 0, error('N must be odd'); end
if N < 3,         error('N must be at least 3'); end

h    = T / (N - 1);
w    = zeros(1, N);
w(1) = h/3;
w(N) = h/3;

for i = 2:N-1
    if mod(i-1, 2) == 1   % i-1 odd => i even => coefficient 4
        w(i) = 4*h/3;
    else                   % i-1 even => i odd => coefficient 2
        w(i) = 2*h/3;
    end
end
end

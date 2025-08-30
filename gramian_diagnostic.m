function gramian_diagnostic()
% Comprehensive diagnostic of Gramian computation
clc; clear;

fprintf('=== GRAMIAN COMPUTATION DIAGNOSTIC ===\n');

% Use exact paper parameters
n = 2; m = 1; T = 2*pi; N = 101;
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

% Test different versions of K(t)
K_versions = {
    @(t) 14.958 * [1 + 0.2*cos(t); 0.5*sin(t)],  % As written in paper
    @(t) 1.4958 * [1 + 0.2*cos(t); 0.5*sin(t)],  % Without one zero
    @(t) 0.14958 * [1 + 0.2*cos(t); 0.5*sin(t)], % Decimal shift
    @(t) [1 + 0.2*cos(t); 0.5*sin(t)] / 14.958,   % Inverse scaling
    @(t) [1 + 0.2*cos(t); 0.5*sin(t)] * 0.1,      % Small scaling
};

version_names = {
    'Paper exact (14.958)', 'Reduced 1 (1.4958)', 'Reduced 2 (0.14958)', ...
    'Inverse (1/14.958)', 'Small scale (0.1)'
};

fprintf('Testing different K(t) interpretations:\n');
for v = 1:length(K_versions)
    fprintf('\n--- %s ---\n', version_names{v});
    K_func = K_versions{v};
    
    % Test K(t) at a few points
    fprintf('K(0) = [%.4f; %.4f]\n', K_func(0));
    fprintf('K(π/2) = [%.4f; %.4f]\n', K_func(pi/2));
    
    try
        % Compute Gramian
        W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
        
        % Analyze
        sigma_vals = svd(W);
        sigma_min = min(sigma_vals);
        sigma_max = max(sigma_vals);
        kappa = sigma_max / sigma_min;
        
        fprintf('σ_min(W) = %.6e\n', sigma_min);
        fprintf('σ_max(W) = %.6e\n', sigma_max);
        fprintf('κ(W)     = %.6e\n', kappa);
        
        % Check if it matches paper values
        target_sigma = 1.250000e-02;
        target_kappa = 1.059859e+02;
        
        sigma_match = abs(sigma_min - target_sigma) / target_sigma < 0.01;
        kappa_match = abs(kappa - target_kappa) / target_kappa < 0.01;
        
        if sigma_match && kappa_match
            fprintf('✓ EXCELLENT MATCH TO PAPER!\n');
        elseif sigma_match
            fprintf('✓ σ_min matches paper\n');
        elseif kappa_match
            fprintf('✓ κ matches paper\n');
        else
            fprintf('✗ Does not match paper\n');
        end
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
    end
end

% Diagnostic 2: Check the vectorized system matrices
fprintf('\n=== VECTORIZED SYSTEM DIAGNOSTIC ===\n');
K_test = @(t) 0.1 * [1 + 0.2*cos(t); 0.5*sin(t)];

% Check dimensions at t=0
t_test = 0;
A_val = A_func(t_test);
B_val = B_func(t_test);
K_val = K_test(t_test);

fprintf('At t=0:\n');
fprintf('A(0) = \n'); disp(A_val);
fprintf('B(0) = \n'); disp(B_val);
fprintf('K(0) = \n'); disp(K_val);

% Form vectorized system matrices
I_n = eye(n);
A_vec = kron(I_n, A_val) + kron(B_val.', I_n);
K_vec = kron(I_n, K_val);

fprintf('Vectorized system at t=0:\n');
fprintf('A_vec size: %dx%d\n', size(A_vec));
fprintf('K_vec size: %dx%d\n', size(K_vec));
fprintf('A_vec = \n'); disp(A_vec);
fprintf('K_vec = \n'); disp(K_vec);

% Check eigenvalues of A_vec
eig_A = eig(A_vec);
fprintf('Eigenvalues of A_vec: '); disp(eig_A.');

% Diagnostic 3: Simple numerical test
fprintf('\n=== SIMPLE NUMERICAL TEST ===\n');
fprintf('Testing with constant, simple matrices...\n');

% Very simple test case
A_simple = @(t) [0, 1; -1, 0];
B_simple = @(t) zeros(2, 2);
K_simple = @(t) [0.1; 0.1];

try
    W_simple = compute_periodic_gramian_block(A_simple, B_simple, K_simple, T, 21);
    sigma_min_simple = min(svd(W_simple));
    fprintf('Simple case σ_min = %.6e\n', sigma_min_simple);
    
    % Analytical check: for this simple case, we can compute approximately
    % The fundamental matrix for [0,1;-1,0] over [0,2π] should be identity
    % So roughly W ≈ ∫₀^(2π) K K^T dt = 2π * [0.01, 0.01; 0.01, 0.01]
    W_approx = 2*pi * 0.01 * ones(4, 4);  % This is wrong but gives idea
    fprintf('Expected order of magnitude: ~%.2e\n', 2*pi*0.01);
    
catch ME
    fprintf('Simple test failed: %s\n', ME.message);
end

% Diagnostic 4: Check quadrature accuracy
fprintf('\n=== QUADRATURE DIAGNOSTIC ===\n');
fprintf('Testing quadrature convergence...\n');

K_func = @(t) 0.1 * [1 + 0.2*cos(t); 0.5*sin(t)];
N_values = [11, 21, 41, 81, 101];

fprintf('N    σ_min(W)\n');
for N_test = N_values
    W_test = compute_periodic_gramian_block(A_func, B_func, K_func, T, N_test);
    sigma_min_test = min(svd(W_test));
    fprintf('%3d  %.6e\n', N_test, sigma_min_test);
end

fprintf('\n=== DIAGNOSTIC COMPLETE ===\n');
end

% Support function (same as before)
function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
K0 = K_func(0);
[n, m] = size(K0);

if mod(N,2)==0, error('N must be odd for composite Simpson rule'); end
tau = linspace(0, T, N);
w = simpson_weights(N, T);

W = zeros(n^2, n^2);

for i = 1:N
    Ki = K_func(tau(i));
    M_i = zeros(n^2, m*n);

    for k = 1:m
        zcol = Ki(:, k);
        for j = 1:n
            ej = zeros(n,1); ej(j)=1;
            Z0 = zcol * (ej.');

            if abs(tau(i) - T) < 1e-10
                Z_final = Z0;
            else
                sylv_ode = @(t, Z) sylvester_rhs(t, Z, A_func, B_func);
                opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
                [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
                Z_final = reshape(Z_sol(end, :), n, n);
            end

            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end

    W = W + w(i) * (M_i * M_i');
end
end

function dZ = sylvester_rhs(t, Z_vec, A_func, B_func)
n = round(sqrt(length(Z_vec)));
Z = reshape(Z_vec, n, n);
dZ = A_func(t) * Z + Z * B_func(t);
dZ = dZ(:);
end

function w = simpson_weights(N, T)
if mod(N,2)==0, error('N must be odd for composite Simpson rule'); end
h = T/(N-1);
w = zeros(1,N);
w(1) = h/3; w(N) = h/3;
for i = 2:N-1
    if mod(i-1,2)==0, w(i) = 4*h/3; else, w(i) = 2*h/3; end
end
end

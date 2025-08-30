function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%COMPUTE_PERIODIC_GRAMIAN_BLOCK Block-wise computation of reachability Gramian
%   Computes the reachability Gramian for periodic Sylvester matrix systems
%   using a structure-exploiting block propagation method that avoids forming
%   large Kronecker matrices.
%
%   INPUTS:
%     A_func - Function handle for A(t) (n x n matrix)
%     B_func - Function handle for B(t) (n x n matrix)  
%     K_func - Function handle for K(t) (n x m matrix)
%     T      - Period of the system
%     N      - Number of quadrature nodes (must be odd for Simpson's rule)
%
%   OUTPUT:
%     W      - Reachability Gramian (n^2 x n^2 matrix)
%
%   ALGORITHM:
%     The algorithm exploits the Sylvester structure by solving n*m individual
%     n x n matrix ODEs instead of one large n^2 x n^2 system. This reduces
%     complexity from O(N*n^6) to O(N*n^3*m).
%
%   REFERENCE:
%     M. S. V. D. Sudarsan, "Controllability and Efficient Gramian Computation 
%     for Periodic Sylvester Matrix Systems", Applied Mathematics Letters (2025)

% Input validation
if nargin < 5
    error('Five input arguments required: A_func, B_func, K_func, T, N');
end

if mod(N, 2) == 0
    error('N must be odd for composite Simpson''s rule');
end

if N < 3
    error('N must be at least 3 for Simpson''s rule');
end

% Get system dimensions
try
    K0 = K_func(0);
    [n, m] = size(K0);
    
    A0 = A_func(0);
    B0 = B_func(0);
    
    if size(A0, 1) ~= n || size(A0, 2) ~= n
        error('A(t) must be n x n where n is the number of rows in K(t)');
    end
    
    if size(B0, 1) ~= n || size(B0, 2) ~= n
        error('B(t) must be n x n where n is the number of rows in K(t)');
    end
    
catch ME
    error('Error evaluating system matrices at t=0: %s', ME.message);
end

% Setup quadrature (composite Simpson's rule)
tau = linspace(0, T, N);
w = simpson_weights(N, T);

% Initialize Gramian
W = zeros(n^2, n^2);

% Progress indicator for large computations
if n*m*N > 1000
    fprintf('Computing Gramian: %d quadrature nodes, %d ODEs per node...\n', N, n*m);
end

% Main loop over quadrature nodes
for i = 1:N
    % Evaluate K at current quadrature node
    Ki = K_func(tau(i));
    
    % Initialize propagation matrix for this quadrature node
    M_i = zeros(n^2, m*n);
    
    % Loop over input columns
    for k = 1:m
        % Extract k-th column of K(tau_i)
        k_col = Ki(:, k); % n x 1 vector
        
        % Loop over basis columns (standard basis for R^n)
        for j = 1:n
            % Create j-th standard basis vector
            e_j = zeros(n, 1);
            e_j(j) = 1;
            
            % Initial condition: Z0 = k_col * e_j^T
            % This creates an n x n matrix with k_col in the j-th column
            Z0 = k_col * e_j.';
            
            % Handle the endpoint case (tau_i = T)
            if abs(tau(i) - T) < 1e-10
                Z_final = Z0;
            else
                % Solve the Sylvester matrix ODE: dZ/dt = A(t)*Z + Z*B(t)
                % from t = tau_i to t = T with initial condition Z(tau_i) = Z0
                try
                    sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func, n);
                    
                    % ODE solver options for high accuracy
                    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, ...
                                  'MaxStep', (T-tau(i))/10);
                    
                    % Solve ODE (vectorized form)
                    [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
                    
                    % Reshape solution back to matrix form
                    Z_final = reshape(Z_sol(end, :), n, n);
                    
                catch ME
                    error('ODE solver failed at node %d, input %d, basis %d: %s', ...
                          i, k, j, ME.message);
                end
            end
            
            % Store the vectorized result in the appropriate column of M_i
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end
    
    % Accumulate contribution to Gramian: W += w_i * M_i * M_i^T
    W = W + w(i) * (M_i * M_i');
    
    % Progress update for large computations
    if n*m*N > 1000 && mod(i, max(1, floor(N/10))) == 0
        fprintf('  Progress: %d/%d nodes (%.1f%%)\n', i, N, 100*i/N);
    end
end

% Ensure Gramian is exactly symmetric (numerical cleanup)
W = (W + W') / 2;

% Verify positive semidefiniteness
min_eig = min(eig(W));
if min_eig < -1e-12
    warning('Gramian has significantly negative eigenvalue: %.2e', min_eig);
end

if n*m*N > 1000
    fprintf('Gramian computation complete.\n');
end

end

function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func, n)
%SYLVESTER_RHS Right-hand side for Sylvester matrix ODE
%   Computes dZ/dt = A(t)*Z + Z*B(t) in vectorized form
%   
%   INPUT:
%     t      - Current time
%     Z_vec  - Current state (vectorized n x n matrix)
%     A_func - Function handle for A(t)
%     B_func - Function handle for B(t)
%     n      - Matrix dimension
%
%   OUTPUT:
%     dZ_vec - Time derivative (vectorized form)

% Reshape vectorized state back to matrix form
Z = reshape(Z_vec, n, n);

% Evaluate system matrices at current time
At = A_func(t);
Bt = B_func(t);

% Compute Sylvester equation: dZ/dt = A(t)*Z + Z*B(t)
dZ = At * Z + Z * Bt;

% Return in vectorized form
dZ_vec = dZ(:);

end

function w = simpson_weights(N, T)
%SIMPSON_WEIGHTS Composite Simpson's rule quadrature weights
%   Computes weights for composite Simpson's rule integration over [0,T]
%
%   INPUT:
%     N - Number of nodes (must be odd)
%     T - Integration interval length
%
%   OUTPUT:
%     w - Quadrature weights (1 x N vector)

if mod(N, 2) == 0
    error('N must be odd for composite Simpson''s rule');
end

if N < 3
    error('N must be at least 3 for Simpson''s rule');
end

% Step size
h = T / (N - 1);

% Initialize weights
w = zeros(1, N);

% Endpoint weights
w(1) = h/3;
w(N) = h/3;

% Interior weights: 4h/3 for odd indices, 2h/3 for even indices
for i = 2:N-1
    if mod(i-1, 2) == 1  % i-1 is odd (i is even)
        w(i) = 4*h/3;
    else                 % i-1 is even (i is odd)
        w(i) = 2*h/3;
    end
end

end

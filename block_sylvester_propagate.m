function M_i = block_sylvester_propagate(A_func, B_func, K_at_tau, t0, tf, varargin)
% BLOCK_SYLVESTER_PROPAGATE
% Build the block M_i \in R^{n^2 x (m n)} at time t0 for the Gramian integrand,
% by propagating Z via dZ/dt = A(t) Z + Z B(t) from t0 to tf with
% initial conditions Z0 = K(:,k) e_j^T for k=1..m and j=1..n.
%
% Usage:
%   M_i = block_sylvester_propagate(A_func, B_func, K_at_tau, t0, tf)
%   M_i = block_sylvester_propagate(..., 'RelTol',1e-9, 'AbsTol',1e-12, 'Solver','ode45')
%
% Inputs:
%   A_func    - handle: A(t) -> n x n
%   B_func    - handle: B(t) -> n x n
%   K_at_tau  - n x m matrix, K(t0)
%   t0, tf    - start and end times (with t0 <= tf)
%
% Optional name-value:
%   'RelTol'  - ODE relative tolerance (default 1e-9)
%   'AbsTol'  - ODE absolute tolerance (default 1e-12)
%   'Solver'  - ODE solver name (default 'ode45', alternative 'ode15s')
%
% Output:
%   M_i       - n^2 x (m n) matrix; columns are vec(Z_final^(k,j)) for each (k,j)
%
% This function implements the block-wise propagation step consistent with
% the repository’s Gramian computation and the paper’s algorithmic description.
%
% Author:  M. S. V. D. Sudarsan
% Email:   msvdsudarsan@gmail.com
% Date:    August 2025

    % Dimensions
    [n, m] = size(K_at_tau);

    % Defaults
    p = inputParser;
    addParameter(p, 'RelTol', 1e-9, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'AbsTol', 1e-12, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'Solver', 'ode45', @(s)ischar(s) || isstring(s));
    parse(p, varargin{:});
    reltol = p.Results.RelTol;
    abstol = p.Results.AbsTol;
    solver = char(p.Results.Solver);

    % Preallocate output
    M_i = zeros(n^2, m*n);

    % Choose solver
    use_ode15s = strcmpi(solver,'ode15s');

    % Loop over input columns and basis columns
    for k = 1:m
        Kk = K_at_tau(:, k); % n x 1

        for j = 1:n
            ej = zeros(n,1); ej(j) = 1;

            % Initial matrix: only j-th column equals K(:,k)
            Z0 = Kk * (ej.');

            % Vectorize and integrate
            z0_vec = Z0(:);
            rhs = @(t,z) sylvester_rhs_vec(t, z, A_func, B_func, n);
            opts = odeset('RelTol', reltol, 'AbsTol', abstol);

            if use_ode15s
                [~, zsol] = ode15s(rhs, [t0 tf], z0_vec, opts);
            else
                [~, zsol] = ode45(rhs, [t0 tf], z0_vec, opts);
            end

            Z_final = reshape(zsol(end, :).', n, n);
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end
end

function dz_dt = sylvester_rhs_vec(t, z_vec, A_func, B_func, n)
    Z = reshape(z_vec, n, n);
    dZ = A_func(t) * Z + Z * B_func(t);
    dz_dt = dZ(:);
end

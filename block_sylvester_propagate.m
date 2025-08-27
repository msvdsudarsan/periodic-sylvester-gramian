function M_i = block_sylvester_propagate(A_func, B_func, K_at_tau, t0, tf, varargin)
% BLOCK_SYLVESTER_PROPAGATE
% Build M_i \in R^{n^2 x (m n)} at time t0 by propagating
% dZ/dt = A(t) Z + Z B(t) from t0 to tf with Z0 = K(:,k) e_j^T.

    [n, m] = size(K_at_tau);

    p = inputParser;
    addParameter(p, 'RelTol', 1e-9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'AbsTol', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'Solver', 'ode45', @(s)ischar(s)||isstring(s));
    parse(p, varargin{:});
    reltol = p.Results.RelTol;
    abstol = p.Results.AbsTol;
    solver = char(p.Results.Solver);

    M_i = zeros(n^2, m*n);
    use_ode15s = strcmpi(solver,'ode15s');

    for k = 1:m
        Kk = K_at_tau(:, k);
        for j = 1:n
            ej = zeros(n,1); ej(j)=1;
            Z0 = Kk * (ej.');
            z0_vec = Z0(:);
            rhs = @(t,z) sylvester_rhs_vec(t, z, A_func, B_func, n);
            opts = odeset('RelTol',reltol,'AbsTol',abstol);
            if use_ode15s
                [~, zsol] = ode15s(rhs, [t0 tf], z0_vec, opts);
            else
                [~, zsol] = ode45(rhs, [t0 tf], z0_vec, opts);
            end
            Z_final = reshape(zsol(end, :).', n, n);
            M_i(:, (k-1)*n + j) = Z_final(:);
        end
    end
end

function dz_dt = sylvester_rhs_vec(t, z_vec, A_func, B_func, n)
    Z = reshape(z_vec, n, n);
    dZ = A_func(t) * Z + Z * B_func(t);
    dz_dt = dZ(:);
end

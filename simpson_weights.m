function w = simpson_weights(N, T)
% SIMPSON_WEIGHTS Generate weights for composite Simpson's rule
%
% Inputs:
%   N - Number of points (must be odd)
%   T - Interval length [0, T]
%
% Output:
%   w - Weight vector (N x 1)

    if mod(N, 2) == 0
        error('N must be odd for Simpson''s rule');
    end
    
    h = T / (N - 1);
    w = zeros(N, 1);
    
    % Simpson's 1/3 rule weights
    w(1) = h/3;
    w(N) = h/3;
    
    for i = 2:N-1
        if mod(i, 2) == 0
            w(i) = 4*h/3;  % Even indices get weight 4h/3
        else
            w(i) = 2*h/3;  % Odd indices get weight 2h/3
        end
    end
end

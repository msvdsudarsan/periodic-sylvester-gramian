function w = simpson_weights(N, T)
% Composite Simpson weights (N must be odd)
if mod(N,2)==0
    error('N must be odd for composite Simpson rule');
end
h = T/(N-1);
w = zeros(1,N);
w(1) = h/3; w(N) = h/3;
for i = 2:N-1
    if mod(i-1,2)==0, w(i) = 4*h/3; else, w(i) = 2*h/3; end
end
end

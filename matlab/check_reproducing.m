% Check the reproduction of parabolas in R^n, n>1, on non-uniform grids

t = [-2,0,1,4];

% f = @(t) [t.^2;1+t+t.^2]; % n = 2
% f = @(t) [t.^2;1+t+t.^2;(1-3*t).^2+5]; % n = 3
f = @(t) [t.^2;1+t+t.^2;(1-3*t).^2+5;-t-4*t.^2]; % n = 4
lvl = 1;

[ f1 ] = subdivision_scheme( f(t),lvl,2 );
f1 - f([0.25,0.75]) % must be zero up to machine precision
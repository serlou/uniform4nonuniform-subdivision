% Just run this

close all;

r = 1;
N = 10;
T = 2*pi;
h = T/(N-1);

u = 0:h:T;
u = u(1:end-1);
v = 0;

[x,y,z] = trefoil(u, v, r);

f = [x;y;z];
example_aux(f);
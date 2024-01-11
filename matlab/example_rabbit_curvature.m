% Just run this

close all;
load('rabbit_data.mat');

f = [f,f];
f1 = subdivision_scheme(f,12,2);

% compute the curvature
[L,R,k] = curvature(f1');
figure;
plot(L(1:end/2),1./R(1:end/2),'.k','MarkerSize',8);
xlabel('Cumulative arc length')
ylabel('Curvature')
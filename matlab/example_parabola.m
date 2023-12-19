% Just run this

close all;
n_steps = 5;
fun = @(x) (x).^2;
% t = 4*(rand(1,15)-0.5);
% t = sort(t);
t = [-1.9381   -1.7893   -1.5751   -1.3313   -1.2075   -0.9235   -0.8321   -0.6420   -0.5104   -0.2734   -0.0412    0.9514    1.6813    1.8065    1.9363];
dt = diff(t);
alpha = dt(1:end-2)./dt(2:end-1);
beta = dt(3:end)./dt(2:end-1);
rho = max([min(alpha)^(-1)-1
max(alpha)-1
min(beta)^(-1)-1
max(beta)-1]);
disp(['Minimum value of \rho reproducing the parabola on this grid: ',num2str(rho)])
f = [t;fun(t)];

figure;
hold on;

show = @(f,opt) plot(f(1,:), f(2,:), opt, 'MarkerSize', 15, ...
    'LineWidth', 1.5);

show(f,'.k');

f0 = subdivision_scheme(f,n_steps,0);
show(f0,':b');

f2 = subdivision_scheme(f,n_steps,2);
show(f2,'-g');

f6 = subdivision_scheme(f,n_steps,6);
show(f6,'--r');

legend('Data','\rho=0','\rho=2','\rho=6');

h = figure;
hold on;
t_final = grid_scheme(t,n_steps);
p_final = fun(t_final);
show_error = @(f,opt) plot(t_final, abs(f(2,:)-p_final), opt, ...
    'MarkerSize', 15, 'LineWidth', 1.5);
show_error(f0,':b');
show_error(f2,'-g');
show_error(f6,'--r');
h.Children(1).YScale = 'log';

title('Approximation error')
legend('\rho=0','\rho=2','\rho=6');

function [ f1 ] = grid_scheme( f,lvl )
% Chaikin's subdivision scheme for data in  R^n
% INPUT: f - matrix of size n x m, where m is the number of points in R^n, n>=1
%        n_steps - refinement steps
% OUTPUT: f1 - matrix of size n x ( 6 + 2^n_steps ( m - 6 ) ), which is the
%              result of n_steps refinement steps of the input data f
if lvl>0
    f1 = zeros(size(f,1),2*size(f,2)-6);
    for i=1:size(f,2)-3
        ff = f(:,i:i+3);
        f1(:,2*i-1) = (3*ff(:,2)+ff(:,3))/4;
        f1(:,2*i) = (ff(:,2)+3*ff(:,3))/4;
    end
    f1 = grid_scheme(f1,lvl-1);
else
    f1 = f;
end
end
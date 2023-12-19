function example_aux(f)
% Auxilary function to generate the examples
% INPUT: f - matrix of size n x m, where m is the number of points in R^n,
%           n = 2 or n = 3, describing a closed curve. f(:,1) must coincide
%           with f(:,end)

n_steps = 5;    % refinement levels

figure;
hold on;
if size(f,1) == 2
    show = @(f,opt) plot(f(1,:), f(2,:), opt, 'MarkerSize', 15, ...
        'LineWidth', 1.5);
    axis image;
    axis off;
else
    show = @(f,opt) plot3(f(1,:), f(2,:), f(3,:), opt, 'MarkerSize', 15,...
        'LineWidth', 1.5);
    view([45,45])
    axis equal;
    axis on;
end

show(f,'.k');

f = [f,f];


f1 = subdivision_scheme(f,n_steps,0);
show(f1,':b');

f1 = subdivision_scheme(f,n_steps,1);
show(f1,'-.m');

f1 = subdivision_scheme(f,n_steps,2);
show(f1,'-g');

f1 = subdivision_scheme(f,n_steps,6);
show(f1,'--r');

legend('Data','\rho=0','\rho=1','\rho=2','\rho=6');
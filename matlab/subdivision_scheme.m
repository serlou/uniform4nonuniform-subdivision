function [ f1 ] = subdivision_scheme(f, n_steps, rho)
% Nonlinear subdivision scheme for data in  R^n
% INPUT: f - matrix of size n x m, where m is the number of points in R^n, n>=2
%        n_steps - refinement steps
%        rho - parameter in [0,+inf). rho = 0 corresponds to linear scheme
% OUTPUT: f1 - matrix of size n x ( 6 + 2^n_steps ( m - 6 ) ), which is the
%              result of n_steps refinement steps of the input data f

if n_steps>0
    f1 = zeros(size(f,1),2*size(f,2)-6); % preallocate memory
    for i=1:size(f,2)-3 
        ff = f(:,i:i+3);    % for each 4 points
        dff = diff(ff,1,2);
        % apply the Lagrange linear rule 'Lambda' but with alpha and beta as functions of data
        f1(:,[2*i-1,2*i]) = Lambda(alpha_beta(rho,A(dff),B(dff)),ff);
    end
    f1 = subdivision_scheme(f1,n_steps-1,rho/2); % recursive call. rho is halved at each step
else
    f1 = f; % if n_steps = 0, return the same data
end
end

function x = A(delta)
    num = det([
        delta(:,1)'*delta(:,2) delta(:,1)'*delta(:,3)
        delta(:,3)'*delta(:,2) delta(:,3)'*delta(:,3)]);
    den = det([
        delta(:,1)'*delta(:,1) delta(:,1)'*delta(:,3)
        delta(:,1)'*delta(:,3) delta(:,3)'*delta(:,3)]);
    if num*den == 0
        x = 1/2;
    else
        x = abs(num/den);
    end        
end

function x = B(delta)
    x = A(delta(:,[3,2,1]));
end

function x = alpha(rho,A,B)
    x = max(1/(1+rho), min(1+rho, 1/(A+sqrt(A*(A+1)*B*(B+1))/(B+1)) ));
end

function x = beta(rho,A,B)
    x = alpha(rho,B,A);
end

function xx = alpha_beta(rho,A,B)
    xx = [alpha(rho,A,B),beta(rho,A,B)];
end

function f1 = Lambda(alpha_beta,f)
    % Lagrange linear rule
    
    alpha = alpha_beta(1);
    beta = alpha_beta(2);
    f1 = zeros(size(f,1),2);
    f1(:,1) = am1(alpha,beta)*f(:,1) + a0(alpha,beta)*f(:,2) + a1(alpha,beta)*f(:,3) + a2(alpha,beta)*f(:,4);
    f1(:,2) = a2(beta,alpha)*f(:,1) + a1(beta,alpha)*f(:,2) + a0(beta,alpha)*f(:,3) + am1(beta,alpha)*f(:,4);
end

function x = am1(alpha,beta)
x = -(3*(4*beta+3))/(64*(alpha^2+alpha)*(alpha+beta+1));
end

function x = a0(alpha,beta)
x = (3*(4*alpha+1)*(4*beta+3))/(64*alpha*(beta+1));
end

function x = a1(alpha,beta)
x = ((4*alpha+1)*(4*beta+3))/(64*beta*(alpha+1));
end

function x = a2(alpha,beta)
x = -(3*(4*alpha+1))/(64*(beta^2+beta)*(alpha+beta+1));
end
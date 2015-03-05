%% Example of linear CG method for a matrix with clustered eigenvalues

clc, clear all, clf

% This is the implementation of the 
% conjugate gradient method (Algorithm 5.1 in
% Nocedal and Wright "Numerical Optimization") 
% for solving linear systems Ax=b, with A symmetric and 
% positive definite matrix, 

A=[];
for i=1:10
    B=rand(10); B=B'*B;
    B=B/max(eig(B))+i*eye(10);
    A=blkdiag(A,B);
end

b=rand(100,1);

n=length(b);
xstar = A\b; % uses the 'exact' solution for computing the error 
% done at each iteration

% choose the number of iterations to be performed
Niter = n;

res = zeros(Niter,1);
err = zeros(Niter,1);
errA = zeros(Niter,1);


x = zeros(n,1); % initial guess
r = A*x-b;

p = -r; % initial direction (steepest descent)
for k = 1:Niter
    alpha = -(r'*p)/(p'*A*p);
    x = x + alpha*p;
    r = A*x-b;
    res(k) = r'*A*r; % computes the A-norm of the residuals
   % err(k) = norm(xstar-x); % computes the errors
    errA(k) = (x-xstar)'*r; % computes the error in the A-norm
    beta = (r'*A*p)/(p'*A*p);
    p = -r+beta*p;  
end
xCG = x;

subplot(3,1,1)
plot(eig(A),zeros(100,1),'o')

subplot(3,1,2)
semilogy(res,'-'); grid on
title('Log Plot of the residual:  $\| r^k\|_A^2$ ', 'interpreter','latex') 

subplot(3,1,3)
semilogy(errA,'-o'); grid on
xlabel('iteration k')
title('Log Plot of the error:  $\log \| x^k-x^*\|_A^2$','interpreter','latex') 

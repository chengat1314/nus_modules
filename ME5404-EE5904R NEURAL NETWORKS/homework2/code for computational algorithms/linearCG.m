%% Linear Conjugate Gradient Method
% this is a simple implementation of the 
% conjugate gradient method used for solving 
% linear systems Ax=b

A=[1 2;2 5]; % must be positive definite symmetric matrix
b=[2;-3];
n=length(b);

x=[0;0]; % initial search
r=A*x-b; 
p=-r; % initial direction (steepest descent)
for k=1:n
    alpha = -(r'*p)/(p'*A*p);
    x=x+alpha*p;
    r=A*x-b;
    beta = (r'*A*p)/(p'*A*p);
    p=-r+beta*p;
    x
end

% comparison with the exact solution 
% found via Gaussian elimination
xstar=A\b 

error = norm(xstar-x)
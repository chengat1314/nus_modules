% This code implements the Fletcher-Reeves method (nonlinear CG)
% for the Rosenbrock function

tol=10^(-6);
n=2;
x=zeros(n,1); % initial search
[y,grad,hy]=rosenbrock(x);
r=grad;
res=norm(r);
p=-r; % initial direction (steepest descent)

while norm(r)>tol
    alpha=linesearch_secant(@rosenbrock,p,x)
    [y,grad,hy]=rosenbrock(x); 
    x=x+alpha*p;
    [y,gradnew,hy]=rosenbrock(x);
    r=gradnew;
    res=[res;norm(r)];
    beta = (norm(gradnew))^2/(norm(grad)^2);
    p=-r+beta*p; 
end
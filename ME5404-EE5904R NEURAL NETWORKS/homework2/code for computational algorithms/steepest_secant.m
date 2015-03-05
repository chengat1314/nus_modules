function [x, xiter, ithist, iflag] = steepest_secant( f, x, tol, maxit)
%
%  function [x, xiter, ithist, iflag] = newton( f, x, tol, maxit )
%
%  newton is the n dimensional Newton's method that solves
%  Min f(x), where f is a n dimensional multivariable function and
%  x is a n dimensional vector
%  g(x) is the gradient of f -- n dimensional vector
%  H(x) is the Hessian of f -- n by n symmetric matrix
%  No line search -- We take the full Newton step in each iteration
%
%  Input parameters:
%    f       name of a matlab function that evaluates
%            f, gradient of f, and the Hessian of f
%    x      Starting solution
%    tol    stopping tolerance (optional. Default tol = 1.e-7)
%            Newton's method stops if  ||g(x)|| < tol
%            where s = -H(x)^{-1} g(x)  is the Newton step.
%    maxit   maximum number of iterations (optional. Default maxit = 100)
%
%
%  Output parameters:
%    x       Approximation to the solution
%    xiter   History of iterates.
%    ithist  array with the iteration history
%            The i-th row of ithist contains  [it, ||f(x^i)||, ||s^i||, ||x^{i+1} - x^i||,
%    ||x^{i+1} - x^i||/||x^i - x^{i-1}||]
%    ifag    return flag
%            iflag =  0  ||f(x)|| <= tolf
%            iflag =  1  iteration terminated because maximum number of
%                        iterations was reached. ||g(x)|| > tol
%

% set tolerances if necessary
if( nargin < 3 ) tol = 1e-7; maxit = 100; end
if( nargin < 4 ) maxit = 100; end

iflag    = 0;
% fx is the value of f
% gx is the gradient of f
% Hx is the Hessian of f
[fx,gx,Hx] = feval(f, x);
gxnorm   = norm(gx,2);
xiter{1} = x;
% Options for MATLAB's eigenvalue solver 'eigs'
options.disp = 0;
options.isreal = 1;
options.issym = 1;
delta = 1e-3;

it       = 1;
snorm    = 100;
gxnorm   = 100;
err      = 100;
n = length(Hx);

while( it < maxit & gxnorm > tol )
    % Use A\b to solve Ax=b
    D = eigs(Hx,1,'SA',options);
    % If Hx is not positive definite, add a diagonal perturbation
    Hx = Hx + max(0,delta-D)*speye(n);
    %s  = -Hx\gx;
    s = -gx;
    alphak = linesearch_secant(f,s,x);
    x = x + alphak*s;
    % xiter is a structure that contains the previous x vectors
    xiter{it+1} = x;
    diff = norm(xiter{it+1} - xiter{it}, 2);
    err = [err; diff];
    ithist(it,:) = [it, fx, gxnorm, snorm, err(it), err(it+1)/err(it)];
    it = it+1;
    [fx,gx,Hx] = feval(f,x);
    gxnorm = norm(gx,2);
    snorm = norm(s,2);
end


% check why the Newton's method truncated and set iflag
if( gxnorm > tol )
    % Newton's method truncated because the maximum number of iterations
    % was reached
    iflag = 1;
    return
else
    % Newton's method truncated because norm(gx,2) <= tol
    % print info for last iteration
    ithist(it,:) = [it, fx, gxnorm, snorm, err(it), err(it)/err(it-1)];
end

end
%%%%%%%%%%%%%%%%%%
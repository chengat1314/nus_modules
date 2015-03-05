function [fx,gx,Hx] = myfunc(x)
% function [fx,gx,Hx] = myfunc(x)
% My test function for Newton's method
% f(x) = (x_1-2)^4 + (x_1-2x_2)^2
% fx, gx, Hx return the function, gradient, and Hessian values at x

fx = (x(1)-2)^4 + (x(1)-2*x(2))^2;
gx = [4*(x(1)-2)^3 + 2*(x(1)-2*x(2)); -4*(x(1)-2*x(2))];
Hx = [12*(x(1)-2)^2+2 -4; -4 8];
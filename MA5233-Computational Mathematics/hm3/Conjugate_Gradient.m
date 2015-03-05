function [x,k]=Conjugate_Gradient(A,b,n_max)
%% author:chengf 20140907@nus
% Solve the linear equations problem by iteration Conjugate_Gradient based
% on the steepest method theory
% % eg: 
% A = [5,1,1;1,5,1;1,1,5];b = [7,7,7]';n_max = 200,
error = 10^(-10);  
[m,n] = size(A);
x = zeros(m,1); % start point
% result = zeros(n_max,m); % result saving
result1 = zeros(1,m);
result2 = zeros(1,m);
%% initial the first direction and first point
x0 = zeros(m,1);
r1 = b - A*x0; d=r1;
result1 = x0;
%% %CG iteration start from here
k = 1;
for i=2:n_max
    alpha = r1'*r1/(d'*A*d);
    result2 = result1 + alpha*d ;
    r2 = r1 - alpha*A*d;
    beta = (r2'*r2)/(r1'*r1);
    d = r2 + beta*d;
    r1 =r2;
    ite_err = r1'*r1;
    result1 = result2;
    if ite_err < error
        sprintf('CG method iteration finished at %i times',i)
        x = result2';
        k = i;
        break
    end
end


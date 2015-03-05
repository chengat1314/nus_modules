function [x,k]=Steepest_decent(A,b,n_max)
%% author:chengf 20140907@nus
% Solve the linear equations problem by iteration SOR method
% % eg:
% clear,clc
% A = [5,1,1;1,5,1;1,1,5];b = [7,7,7]';n_max = 5,
error = 10^(-10);  
[m,n] = size(A);
x = zeros(m,1); % start point
% result = zeros(n_max,m); % result saving
result1 = zeros(1,m);
result2 = zeros(1,m);
%% initial the first direction and first point
x0 = zeros(m,1);
r1 = b - A*x0;  
result1 = x0;
k = 0;
%% %CG iteration start from here  
for i=2:n_max 
    alpha = r1'*r1/(r1'*A*r1);
    result2 = result1 + alpha*r1 ; 
    r1 = b - A*result2; 
    ite_err = r1'*r1;
    result1 = result2;
    if ite_err < error || (r1'*r1 == 0)
        sprintf('CG method iteration finished at %i times',i)
        x  = result2';
        k = i;
        break
    end
end 

function [x,k]=Conjugate_Gradient_forPoissonEquation(A,b,n_max)
%% author:chengf 20140907@nus
 clear,clc
%% initial data
N = 3;
A = two_dim_poisson(10);clear,clc
ax = 0;bx = 1;
ay = 0;by = 1;
hx = (bx - ax)/N;
hy = (by - ay)/N;
f_xy = @(x,y)(sin(5*pi*x) * sin(7*pi*y)); %target function 
b = zeros(N^2,1); k = 1;
for xi = ax + hx:hx:bx
    for yi = ay + hx:hy:by
        b(k) = hx*hy*f_xy(xi,yi); 
        k = k + 1 ;
    end
end
%% index
n_max = 200;
error = 10^(-10);  [m,n] = size(A); x = zeros(m,1); % start point
M = m^(0.5); %dimension
% result = zeros(n_max,m); % result saving
result1 = zeros(1,m);
result2 = zeros(1,m);
%% initial the first direction and first point
x0 = zeros(m,1);
r1 = b - A*x0; d=r1;
result1 = x0;
%% %CG iteration start from here
k = 1;
for i = 1:M
    %format long;
    if( i == 1)
        break
    elseif (i==M)
        break
    else
        break
    end
    %
end
for i=2:n_max
    alpha = r1'*r1/(d'*A*d);
    result2 = result1 + alpha*d; 
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


% for i=2:n_max
%     alpha = r1'*r1/(d'*A*d);
%     result2 = result1 + alpha*d; 
%     r2 = r1 - alpha*A*d;
%     beta = (r2'*r2)/(r1'*r1);
%     d = r2 + beta*d;
%     r1 =r2;
%     ite_err = r1'*r1;
%     result1 = result2;
%     if ite_err < error
%         sprintf('CG method iteration finished at %i times',i)
%         x = result2';
%         k = i;
%         break
%     end
% end


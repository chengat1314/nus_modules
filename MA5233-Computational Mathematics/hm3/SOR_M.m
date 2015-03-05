function x=SOR_M()%(A,b,x,n_max)
%% author:chengf 20140907@nus
% Solve the linear equations problem by iteration SOR method
% eg:
clear,clc
A = [5,1,1;1,5,1;1,1,5],b = [7,7,7]',n_max = 200,
error = 10^(-20);w = 1.5;
[m,n] = size(A);
x = zeros(m,1); % start point
R = zeros(m,n); L = zeros(m,n);   
D = zeros(m,n); U = zeros(m,n);
result = zeros(n_max,m); % result saving  
%% get the iteration matrix first
for i = 1:m
    for j = 1:n
        if i == j 
            D(i,j) = A(i,j);
        elseif i > j 
            L(i,j) = -A(i,j);
        else
            U(i,j) = -A(i,j);
        end
    end
end 
RJ = inv(D) * (L+U) % iteration matrix for jacobi method
w = 2/(1+(1 - max(eig(RJ))^2)^(1/2))  % choose the best w as the iteration parameter
iter_b = w*inv(D - w*L);
R = inv(D - w*L)*((1-w)*D + w*U);  % iteration matrix

%% %SOR iteration start from here
%k = 1;
for i=2:n_max  
    x = R * x + iter_b * b;  
    result(i,:) = x ; 
    ite_err = sum((result(i,:) - result(i-1,:)).^2)/sum((result(i,:)).^2);
    if ite_err < error
        sprintf('SOR method iteration finished at %i times',i)
        x 
        break
    end
end


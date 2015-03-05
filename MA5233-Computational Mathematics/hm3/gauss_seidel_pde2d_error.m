function out = gauss_seidel_pde2d_error(N)
%% Programming the Gauss Siedel Method
% Generating the solution matrix
% clc,clear
% N = 20;
n = N; m= N;
zz=zeros(n+1,m+1); alpha = 0;
f=@(x,y)(x^2 + y);
tolerance = 10^(-6) ; %tolerance
% Assigning the temperature boundary conditions
h = 1/m;
iteration=0;
count=(m-1)*(n-1);
error = zeros(n+1,m+1);
while(count>0)
count=0;
iteration=iteration+1; 
    for i=2:1:n
        for j=2:1:m
            temp=zz(i,j);
            zz(i,j)=(zz(i+1,j)+zz(i-1,j)+zz(i,j+1)+zz(i,j-1) + h^2*f((i-1)*h,(j-1)*h) )/(4+alpha*h^2);
            error(i,j)= abs((zz(i,j)-temp)); 
        end 
    end
    if(sum(sum(error))>tolerance)
        count=1;
    end    
    maxerror(iteration) = max(max(error)); 
end 
out = maxerror(iteration);

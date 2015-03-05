function A = one_dim_poisson(n)
% product one dimensinal poisson matrix 
% n = 4;
A = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = 2;
        elseif i == j-1 || j == i - 1
            A(i,j) = -1;
        end
    end
end



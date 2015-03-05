function A = two_dim_poisson(n)
% product one dimensinal poisson matrix 
%  n = 4;
 A = kron(eye(n),one_dim_poisson(n)) + ...
     kron(one_dim_poisson(n),eye(n));
 

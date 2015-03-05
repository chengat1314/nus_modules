function question3()
%% by least square methods
A=[4,-1,1;2,5,2;1,2,4;4,1,2;2,4,-1;1,1,-3];
b=[8,3,11,9,-5,-9]';
x = inv(A'*A)*A'*b;

%% by SVD methods

[U,S,V] = svd(A) %X = U*S*V'.

x= V*[inv(S(1:3,:)),zeros(3,3)]*U'*b

%% GMRES method
B = [0,0,0,0,1;1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0]
c = [1,0,0,0,0]'

[m,n] = size(B);
x0 = zeros(m,1);
% r = B*x0 - c;
% h = sum(r^2);
GMRES(B,b,x0)
function matrix()
clc,clear;
N =5000;
% a=rand(N,N);
% b=rand(N,N);
% c=zeros(N,N);
% tic;
% for i = 1:N
%     for j = 1:N
%         for k = 1:N
%             c(i,j) = a(i,k) * b(k,j);
%         end
%     end
% end
% toc;
% tic;
% a*b;
% toc;
a=[         0    0.0255    0.3525    0.6669    0.9631
    0.8383    0.3354    0.9153    0.7959    0.8327
    0.3450    0.8712    0.0899    0.8883    0.7010
    0.7346    0.3002    0.0497    0.9082    0.0977
    0.0403    0.0850    0.5588    0.9265    0.0756];
c=zeros(5,5);n=5
for i=1:n
    for j=1:n
    c(i,j) = c(i,j)+a(i,j)*a(j,i);
    end
end
c


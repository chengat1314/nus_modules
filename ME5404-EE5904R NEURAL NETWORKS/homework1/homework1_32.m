function homework1_32()
%homework 1.3 for neural networks
%preceptron
clc,clear,close all;
%% data 
input = [0;1]; %% for COMPLEMENT
output = [1,0];
[m,n]= size(input);
xx = [input,ones(m,1)];
%% preceptron
ita = 1;
w = zeros(n+1,1);
N= 10;
error_his = zeros(N,1);
w_his = zeros(N*m,1);
step = 0;
for k = 1:N
   sumerror = 0;
   for i = 1:m 
      error = output(i)-active_function(xx(i,:)*w); %% rounds each element of t to the nearest number of the specified unit of time greater than or equal to that elemen
       if error~=0
           w = w + ita*error*xx(i,:)';
           sumerror = sumerror + 1 ;
           step = step + 1 ;
           w_his(step) = norm(w);
       end
   end
   error_his(k)=sumerror;
end
w'
error_his

%% plot
figure(1)
plot([1:N*m]',w_his,'.')
xlabel('n');ylabel('norm(w)');
title('n ~ norm(w)');
figure(2)
for i = 1:m 
    if output(i)==0
        plot(input(i,1),0,'r*'); hold on
    else
        plot(input(i,1),0,'gO'); hold on
    end
end
xx = [-0.5:0.01:1.5];
plot(0,xx ,'-r');
xlim([-0.5,1.5])
xlabel('x1');ylabel('x2');
title('x1~x2 Perceptron');
function re = active_function(v)
re = 0;
if v>=0
    re = 1;
end

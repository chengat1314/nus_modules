function  out = forward_euler(f,N,a,b,x0)
%% forward euler method to solve first order ODE
% clc,clear;
% N=30; a=0;b=1;x0 = 1;% sample
% f=@(y,t)(-2*y+1);

k = 1/N; %step length
x = [a:k:b-k];
out = zeros(N,1); out(1) = x0;
%% start iteration 
for i = 2:N
    out(i) = out(i-1) + k*(f(out(i-1),x(i)));
end

%% check the result
% y=0.5*(1+exp(-2.*x));
% plot(x,out,'r*');hold on
% plot(x,y,'y-')

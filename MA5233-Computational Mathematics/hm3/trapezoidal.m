function  out = trapezoidal(f,N,a,b,x0)
%% trapezoidal method to solve first order ODE
% clc,clear;
% N=30; a=0;b=1;x0 = 1;% sample
% f=@(y,t)(-2*y+1);

k = 1/N; %step length
x = [a:k:b-k];
out = zeros(N,1); out(1) = x0;
%% start iteration 
for i = 2:N
    ff = @(y)(y-out(i-1)-(f(out(i-1),x(i-1))+f(y,x(i)))*k/2);
    solu = fzero(ff,0);
    out(i) = solu;
end

%% check the result
% y=0.5*(1+exp(-2.*x));
% plot(x,out,'r*');hold on
% plot(x,y,'y-')

function  out = fourth_runge_kutta(f,N,a,b,x0)
%% 4th order runge kutta method to solve first order ODE
% clc,clear;
% N=30; a=0;b=1;x0 = 1;% sample
% f=@(y,t)(-2*y+1);
% f=@(y,t)(-2*1+t); check variable t
k = 1/N; %step length
x = [a:k:b-k];
out = zeros(N,1); out(1) = x0;
%% start iteration 
for i = 2:N
    k1 = f(out(i-1),x(i-1));
    k2 = f(out(i-1) + k1*k/2,x(i-1)+k/2);
    k3 = f(out(i-1) + k2*k/2,x(i-1)+k/2);
    k4 = f(out(i-1) + k3*k,x(i-1)+k);
    out(i) = out(i-1) + (k/6)*(k1 + 2*k2 + 2*k3 + k4);
end

%% check the result
% y=0.5*(1+exp(-2.*x));
% % y = 1+x.^2/2-2*x;
% plot(x,out,'r*');hold on
% plot(x,y,'y-')

function solve_question_q4()
%% 4th Runge Kutta
clc,clear
N = 50; a=0;b=3;
k=(b-a)/N;
fxy = @(x,y,t)(2*x - y^2 + t*exp(t) - t);
x0=1;  y0=0;
t = [a:k:b-k];
out_y = zeros(N,1); out_y(1) = y0;
out_x = zeros(N,1); out_x(1) = x0;
for i=2:N
    out_y(i) = out_y(i-1) + out_x(i-1)*k;
    
    xk1 = fxy(out_x(i-1),out_y(i-1),t(i-1));
    xk2 = fxy(out_x(i-1) + xk1*k/2,out_y(i-1),t(i-1)+k/2);
    xk3 = fxy(out_x(i-1) + xk2*k/2,out_y(i-1),t(i-1)+k/2);
    xk4 = fxy(out_x(i-1) + xk3*k,out_y(i-1),t(i-1)+k);
    out_x(i) = out_x(i-1) + (k/6)*(xk1 + 2*xk2 + 2*xk3 + xk4); 
end
figure(1)
plot(t,out_x,'g*'); hold on
plot(t,out_y,'r*');


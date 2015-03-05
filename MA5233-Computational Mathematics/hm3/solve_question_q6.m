function solve_question_q6()
%% 4th Runge Kutta
clc,clear
N = 500; a=0;b=4;k=(b-a)/N;
fxy = @(x,y,t)(-4*x+3*y+6*sin(t^2));
gxy = @(x,y,t)(-2.4*x+1.6*y+3.6*cos(t^2));
x0=1;  y0=0;
t = [a:k:b-k];
out_y = zeros(N,1); out_y(1) = y0;
out_x = zeros(N,1); out_x(1) = x0;
for i=2:N
    xk1 = fxy(out_x(i-1),out_y(i-1),t(i-1));
    xk2 = fxy(out_x(i-1) + xk1*k/2,out_y(i-1),t(i-1)+k/2);
    xk3 = fxy(out_x(i-1) + xk2*k/2,out_y(i-1),t(i-1)+k/2);
    xk4 = fxy(out_x(i-1) + xk3*k,out_y(i-1),t(i-1)+k);
    out_x(i) = out_x(i-1) + (k/6)*(xk1 + 2*xk2 + 2*xk3 + xk4); 
    
    yk1 = gxy(out_x(i-1),out_y(i-1),t(i-1));
    yk2 = gxy(out_x(i-1) ,out_y(i-1)+ yk1*k/2,t(i-1)+k/2);
    yk3 = gxy(out_x(i-1) ,out_y(i-1)+ yk2*k/2,t(i-1)+k/2);
    yk4 = gxy(out_x(i-1) ,out_y(i-1)+ yk3*k,t(i-1)+k);
    out_y(i) = out_y(i-1) + (k/6)*(yk1 + 2*yk2 + 2*yk3 + yk4); 
end
figure(1)
plot(t,out_x,'g-'); hold on
% figure(2) 
plot(t,out_y,'r-');



function solve_question_q5()
clc,clear;
f = @(y,t)(t+y);  %-t-1+2*exp(t)
N = 100;a=0;b=1;x0=1; 
x = [a:1/N:b-1/N];
euler = forward_euler(f,N,a,b,x0);
trapez = trapezoidal(f,N,a,b,x0);
fourth_RK = fourth_runge_kutta(f,N,a,b,x0);
exact = [-x-1+2.*exp(x)]';
%% plot numerical result
figure(1)
plot(x,euler,'g*'); hold on
plot(x,trapez,'b+'); hold on
plot(x,fourth_RK,'y-'); hold on
plot(x,exact,'r^'); 
title('numerical solution for different solution');
xlabel('t');
ylabel('y');
legend('euler','trapezoidal','fourth_runge_kutta');

%% plot error
figure(2)
plot(x,log(abs(euler-exact))/log(1/N),'g*'); hold on
plot(x,log(abs(trapez-exact))/log(1/N),'b+'); hold on
plot(x,log(abs(fourth_RK-exact))/log(1/N),'y-');  
xlabel('t');
ylabel('log(error)/log(h)');
legend('euler','trapezoidal','4th runge kutta');
ylim([0 7]);

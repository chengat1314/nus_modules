function res = numerical_sqrt(a,x_initial)
%% athour:chengf@nus
%three numerical methods computing sqrt(a) with different intial values.
% clc,clear
% a = 0.111;
k = 100000;
x0 = x_initial;
ite1=0;ite2=0;ite3=0;
%% newton method base on f=x^2 - a = 0 f'=2x
for i=1:k
    x1 = x0 - (x0^2-a)/(2*x0);
    if abs(x1-x0)<10^(-6)
%         fprintf('the sqrt is %f\n', x1);
%         fprintf('Method1 iter times is %i\n', i);
        break;
    end
    ite1 = i;
    x0=x1;
end
x0 = x_initial;
%% newton method base on f=1 - a/x^2 = 0
for i=1:k
    x1 = x0 - (x0^3-a*x0)/(2*a);
%     x1 = x0 - (1-a/(x0^2))/((2*a)/(x0^3));
    if abs(x1-x0)<10^(-6)
%         fprintf('the sqrt is %f\n', x1);
%         fprintf('Method2 iter times is %i\n', i);
        break;
    end
    ite2 = i;
    x0=x1;
end
x0 = x_initial;
%% third order method
for i=1:k
    x1 = x0*(x0^2 + 3*a)/(3*x0^2 + a);
    if abs(x1-x0)<10^(-6)
%         fprintf('the sqrt is %f\n', x1);
%         fprintf('Method3 iter times is %i\n', i);
        break;
    end
    ite3 = i;
    x0=x1;
end
res = x0;
res = [res,ite1,ite2,ite3];

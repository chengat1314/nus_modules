function integral=trapezoidal()  %(fun,a,b,m,error)
clc
%fun object function; a upper limit; b lower limit; 
%m iteration times; error:precision 
x = 1; error = 1*10^(-10);
%fun = @(t)(t+sin(x*t^2));
fun = @(t)(t^2 + t + 3)
a = 0; b = 1; m = 20;
n=1;h=b-a; T=zeros(1,m+1); t=a;
T(1)=h*(feval(fun,a)+feval(fun,b))/2;
for i=1:m
    h=h/2; n=2*n; s=0;
    for k=1:n/2
        t=a+h*(2*k-1); s=s+feval(fun,t);
    end
    T(i+1)=T(i)/2+h*s;
    T(i+1);
    if abs(T(i+1) - T(i))/abs(T(i))<error
        fprintf('break','%s');i
        break
    end
end

integral=T(i);
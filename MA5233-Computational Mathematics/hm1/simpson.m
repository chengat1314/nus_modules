function y=simpson()  %(fun,a,b,m,error)
clc
%fun object function; a upper limit; b lower limit;
%m iteration times; error:precision
x = 1; error = 1*10^(-6);
%fun = @(t)(t+sin(x*t^2));
fun = @(t)(t^2 + t + 3);
a = 0; b = 1; n = 30;

z1=feval(fun,a)+ feval(fun,b);m=n/2;
h=(b-a)/(2*m); x=a;
z2=0; z3=0; x2=0; x3=0;

for k=2:2:2*m    
    x2=x+k*h; 
    z2= z2+2*feval(fun,x2);    
end

for k=3:2:2*m    
    x3=x+k*h; 
    z3= z3+4*feval(fun,x3);    
end

y=(z1+z2+z3)*h/3;

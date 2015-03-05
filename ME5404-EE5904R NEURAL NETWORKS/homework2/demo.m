function demo()
clc,clear
%% Steep descent
phi = @(x1,x2) 100*(x2 - x1.^2).^2 + (1 - x1).^2;
fx1 = @(x1,x2) -2 + 2*x1 -200*x1*(x2-x1^2);
fx2 = @(x1,x2) 200*(x2-x1^2);

x1 = -1.2;
x2 = 1;
tol = 10;
vold = 10;
k = 1;

while tol > 10^-10
    dx1 = x1 - fx1(x1,x2);
    [x1, vx1] = fminbnd(@(x1) phi(x1,x2), min(x1,dx1), max(x1,dx1));
    dx2 = x2 - fx2(x1,x2);    
    [x2, vx2] = fminbnd(@(x2) phi(x1,x2), min(x2,dx2), max(x2,dx2));
    fprintf('%14.8f %14.8f %14.8f %14.8f\n', x1, vx1, x2, vx2)
    plot(k,x1,'o')
    plot(k,x2,'r*')
    hold on
    tol = abs(vx2-vold);
    vold = vx2;
    k = k+1;
end
legend('x2','x1')
hold off
x1,x2

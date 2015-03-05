function home2_1()
%% E5904R Neural Network
% steepest gradient descent method and newton method
clear,clc;
xy = rand(1,5); % initial point
k = 1; error = 1e-6; ita = 0.001;%learning rate
fxy = @(x,y) 100*(y - x.^2).^2 + (1 - x).^2;
x = xy(1);y=xy(2);
xy(1,3) = fxy(x,y);
%% steepest gradient descent method
df=@(x,y)df(x,y);
while(abs(fxy(x,y))>error & k<50000 )
    pd = df(x,y);
    dx = x - ita*pd(1);
    [x, fx] = fminbnd(@(x) fxy(x,y), min(x,dx), max(x,dx));
    pd = df(x,y);
    dy = y - ita*pd(2);
    [y, fy] = fminbnd(@(y) fxy(x,y), min(y,dy), max(y,dy));
    
    k = k+1;
    xy(k,1:2)=[x,y];
    xy(k,3) = fxy(x,y); %record the values of the f(x,y)
    xy(k,4:5) = df(x,y);
end
xy(end,:)
k
%% plot trajectory of (x,y) in the 2-dimensional space
plotme();
hold on;
x = xy(:,1)';size(x)
y = xy(:,2)';size(y)
plot(x,y,'-*');
title('contour plot of the function and iteration solution');
xlabel('x');
ylabel('y')
%% plot the f(x,y) converges 0
figure(2)
plot(xy(:,3),'-*');
title('the values of f(x,y) with the iteration');
xlabel('step');
ylabel('f(x,y)')
ylim([-0.1,2]);
function fxy = df(x,y) %% differentiation
dfx = 2*(x-1)-400*(y-x^2)*x;
dfy = 200*(y-x^2);
fxy = [dfx,dfy];

function H = Hessian(x,y) %% 
H = [2+1200*x-400*y,-400*x;-400*x,200];

function plotme()
x = linspace(-0.5,2);
y = linspace(-0.5,2);
[X,Y] = meshgrid(x,y);
Z = (1-X).^2+100*(Y-X.^2).^2;
figure(1)
contour(X,Y,Z,50) %,'ShowText','on'



clear all, clc, clf

[X,Y]=meshgrid(-2:0.01:2, -1:0.01:3);

Z=100*(Y-X.^2).^2+(1-X).^2;

v=logspace(0.5,3,30);
contour(X,Y,Z,v);

hold on

plot(1,1,'ro')

%%%
% Steepest descent with computation of step length using the secant method

[x,xiter,ithist,iflag] = steepest_secant(@rosenbrock,[-1.2 1]',1e-6,6000);
fprintf('The solution x computer after %g iterations is', length(xiter))
x

% %%% 
% % This is the steepest descent with symbolic computation of the step
% % length
% 
% [x,xiter,ithist,iflag] = steepest_symbolic(@rosenbrock,[2 1]',1e-6,1000);
% fprintf('The solution x computer after %g iterations is', length(xiter))
% x

for i=2:length(ithist)
    plot([xiter{i}(1) xiter{i-1}(1)], [xiter{i}(2) xiter{i-1}(2)],'o-');
end

hold off





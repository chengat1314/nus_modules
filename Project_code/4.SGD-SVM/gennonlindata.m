function [x, y] = gennonlindata(N, L0, x0, s0, y0);
% function [x, y] = gennonlindata(N, L0, x0, s0, y0);
%
% generate mixture-of-gaussians nonlinear data
%
%     N: dimension of x vectors
%     L0: gaussians' number of points (sums to L)   [C,1] with C = #centers
%     x0: gaussians' centers                        [C,N]
%     s0: gaussians' standard deviations            [C,1]
%     y0: gaussians' labels                         [C,1]
%
%     x: [L,N]
%     y: [L,1] , values in y0, L0 of them each

[C,N] = size(x0);
L = sum(L0);
x = [];
y = [];

for class = (1:C)
    o = ones(L0(class),1);
    x = [x; s0(class)*randn(L0(class),N)+o*x0(class,:)];
    y = [y; y0(class)*o];
end
r = randperm(L);
x = x(r,:);
y = y(r);

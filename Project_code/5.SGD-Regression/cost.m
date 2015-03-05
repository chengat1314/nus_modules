function CostVal = cost(X,y,theta)

% This calculates the cost value. The cost value should always be
% decreasing with each iteration.

% Initialize values
m = length(y); % number of training examples

CostVal = (1/(2*m))*sum((X*theta - y).^2);
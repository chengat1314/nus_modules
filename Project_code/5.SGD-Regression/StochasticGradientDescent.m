function [theta,thetaRecord, CostHistory] = StochasticGradientDescent(X, theta, y, alpha,lambda, numIters)
% StochasticGradientDescent is used to learn the parameters theta in order to fit a
% straight line to the points.

% Initialize values
m = length(y); % number of training examples
CostHistory = zeros(numIters*m, 1);
thetaLen = length(theta);
tempVal = theta; % Just a temporary variable to store theta values.
thetaRecord = zeros(numIters*m,thetaLen);
for iter=1:numIters
    for j = 1:m
        temp = (X(j,:)*theta - y(j));

        for i=1:thetaLen
            tempVal(i,1) = sum(temp.*X(j,i));
        end

        theta(1) =  theta(1)  - (alpha)*tempVal(1);
        theta(2) =  theta(2)*(1-alpha*lambda/m) - (alpha)*tempVal(2);
        CostHistory(iter,1) = cost(X,y,theta);
        thetaRecord(iter,:) = theta;
    end     
end
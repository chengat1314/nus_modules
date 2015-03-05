function [theta,thetaRecord, CostHistory] = gradientDescent(X, theta, y, alpha,lambda, numIters)
% Gradient Descent is used to learn the parameters theta in order to fit a
% straight line to the points.

% Initialize values
m = length(y); % number of training examples
CostHistory = zeros(numIters, 1);
thetaLen = length(theta);
tempVal = theta; % Just a temporary variable to store theta values.
thetaRecord = zeros(numIters,thetaLen);
for iter=1:numIters
    temp = (X*theta - y);
    
    for i=1:thetaLen
        tempVal(i,1) = sum(temp.*X(:,i));
    end
    
%     theta = theta*(1-alpha*lambda/m) - (alpha/m)*tempVal;
    theta(1) =  theta(1)  - (alpha/m)*tempVal(1);
    theta(2) =  theta(2)*(1-alpha*lambda/m) - (alpha/m)*tempVal(2);

    thetaRecord(iter,:) = theta;
    CostHistory(iter,1) = cost(X,y,theta);
 
end
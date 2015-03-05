function homework1_3()
%% Monte Carlo Method
clc,clear
s0= 50; %- initial stock price;
r = 0.05; % risk-free interest rate;
sigma=0.2; %- volatility; 
T=1; %- terminal time; 
K=55; %- strike price;
n1 = 10000;% n - sample size.
m = 50;  % time slots
t=[0: 1/m: 1]*T;  % m monitoring dates 
s = zeros(m,1); % the stock price records
ss= zeros(m,1); % the stock price records
s(1) = s0; 
%% Control Variate Method  
tbar=sum(t(2:m+1))/m;  
mubar=log(s0)+(r-sigma^2/2)*tbar;
sigbar = sqrt(sum((2*m-2*[1:m]+1).*t(2:m+1)))*sigma/m;
theta = (mubar-log(K))/sigbar;
p=(exp(mubar+sigbar^2/2)*normcdf(sigbar+theta)-K*normcdf(theta));
n=2*n1;
for i=1:n    % loops to generate n sample paths and payoffs
  for k=2:m+1   % generate the k-th path and payoff  
    z=randn;  % generate a sample from N(0,1)
    s(k) = s(k-1)*exp((r-0.5*sigma^2)*(t(k)-t(k-1))+sigma*sqrt(t(k)-t(k-1))*z);
  end
  x(i)=max(max(s(2:m+1))-K,0); % payoff for the k-th path
  y(i)=max(geomean(s(2:m+1))-K,0) - p; % control variate
end
xbar = sum(x)/n; ybar = sum(y)/n;
b=sum((x-xbar).*(y-ybar))/sum((y-ybar).^2);% optimal b
% b=0;
H = x-b*y;
v=sum(H)/n; 
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error
fprintf('the result of Control Variate Method:\n v=%8.4e, se=%8.4e\n',v, se);




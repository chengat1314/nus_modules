% Estimate the price of an Asian call option on a stock using the geometric 
% average as the control variate
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_asian_cv(s0, r, sigma, T, m, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; m - number of monitor dates;
% K - strike price;
% n - sample size.

t=[0: 1/m: 1]*T;  % m monitoring dates 
s(1) = s0;

% the expected value of the control variate
tbar=sum(t(2:m+1))/m;  mubar=log(s0)+(r-sigma^2/2)*tbar;
sigbar = sqrt(sum((2*m-2*[1:m]+1).*t(2:m+1)))*sigma/m;
theta = (log(K)-mubar)/sigbar;
p=exp(-r*T)*(exp(mubar+sigbar^2/2)*normcdf(sigbar-theta)-K*normcdf(-theta));

for i=1:n    % loops to generate n sample paths and payoffs
  for k=2:m+1   % generate the k-th path and payoff  
    z=randn;  % generate a sample from N(0,1)
    s(k) = s(k-1)*exp((r-0.5*sigma^2)*(t(k)-t(k-1))+sigma*sqrt(t(k)-t(k-1))*z);
  end
  x(i)=exp(-r*T)*max(sum(s(2:m+1))/m -K,0); % payoff for the k-th path
  y(i)=exp(-r*T)*max(geomean(s(2:m+1))-K,0) - p; % control variate
end

b=0;
%b=1;

% optimal b
xbar = sum(x)/n; ybar = sum(y)/n;
%b=sum((x-xbar).*(y-ybar))/sum((y-ybar).^2);

H = x-b*y;

v=sum(H)/n; 
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('\n v=%8.4e, se=%8.4e\n',v, se);

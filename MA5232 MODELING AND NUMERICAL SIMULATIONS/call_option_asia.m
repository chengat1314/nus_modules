% Example 2: Estimate the price of an Asian call option on a stock
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = call_option_asia(s0, r, sigma, T, m, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; m - number of monitor dates;
% K - strike price;
% n - sample size.

t=[0: 1/m: 1]*T;  % m monitoring dates 
s(1) = s0;

for k=1:n    % loops to generate n sample paths and payoffs
  for i=2:m+1   % generate the k-th path and payoff  
    z=randn;  % generate a sample from N(0,1)
    s(i) = s(i-1)*exp((r-0.5*sigma^2)*(t(i)-t(i-1))+sigma*sqrt(t(i)-t(i-1))*z);
  end
  x(k)=exp(-r*T)*max(sum(s(2:m+1))/m -K,0); % payoff for the k-th path
end

v=sum(x)/n; 
se = sqrt((sum(x.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf(' \n v=%8.4e, se=%8.4e\n',v, se);



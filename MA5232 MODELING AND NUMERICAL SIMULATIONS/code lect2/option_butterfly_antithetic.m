% Estimate the price of a butterfly spread option on a stock
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_butterfly_antithetic(s0, r, sigma, T, K1, k2, k3, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.

z = randn(n,1);  % n samples from N(0,1)

s = s0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*z); % stock price at T
x=exp(-r*T)*(max(s-K1,0)+max(s-k3,0)-2*max(s-k2,0));  % payoff of the option at T

s = s0*exp((r-0.5*sigma^2)*T - sigma*sqrt(T)*z); % antithetic sample
y=exp(-r*T)*(max(s-K1,0)+max(s-k3,0)-2*max(s-k2,0)); 

v = sum(x+y)/(2*n);  % MC estimate of the value of the option
se = sqrt((sum(((x+y)/2).^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('\n v=%8.4e, se=%8.4e\n',v, se);
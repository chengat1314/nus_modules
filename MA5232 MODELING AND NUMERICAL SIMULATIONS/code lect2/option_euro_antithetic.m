% Estimate the price of a European call option using antithetic sampling
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_euro_antithetic(s0, r, sigma, T, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.

z = randn(n,1);  % n samples from N(0,1)

s = s0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*z); % stock price at T
x=exp(-r*T)*max(s-K,0);  % payoff of the option at T

s = s0*exp((r-0.5*sigma^2)*T - sigma*sqrt(T)*z); % antithetic sample
y=exp(-r*T)*max(s-K,0);  % payoff of the option at T

v = sum(x+y)/(2*n);  % MC estimate of the value of the option
se = sqrt((sum(((x+y)/2).^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('\n v=%8.4e, se=%8.4e\n',v, se);
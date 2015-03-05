% Estimate the price of a European call option on a stock using the discounted
%    terminal stock price as the Control Variate
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_euro_cv(s0, r, sigma, T, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.


z = randn(n,1);  % n samples from N(0,1)
s = s0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*z); % stock price at T
x=exp(-r*T)*max(s-K,0);  % payoff of the option at T
y=exp(-r*T)*s-s0; % control Variate

b=0;
b=1;

% optimal b
%xbar = sum(x)/n; ybar = sum(y)/n;
%b=sum((x-xbar).*(y-ybar))/sum((y-ybar).^2);


H=x-b*y;

v = sum(H)/n;  % MC estimate of the value of the option
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('\n v=%8.4e, se=%8.4e\n',v, se);
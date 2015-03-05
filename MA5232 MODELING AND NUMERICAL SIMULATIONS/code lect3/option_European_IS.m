% Estimate the price of an European call option on a stock using improtance sampling
% The stock price is modeled by a geometric Brownian motion.

function [v, se, vexact] = option_European_IS(s0, r, sigma, T, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.


f=@(theta)fun(theta, s0, r, sigma, T, K); % define a function handle f
theta = fzero(f, 10); % find the zero of f; or using fsolve(f, 10)

% theta =0; % plain MC 

y=theta+randn(n,1); % iid samples from the tilted distribution

H=max(s0*exp(-0.5*sigma^2*T+sigma*sqrt(T)*y)-exp(-r*T)*K, 0); %discouted payoff
H=H.*exp(-theta*y+theta^2/2); % discounted payoff x the likelihood ratio

v = sum(H)/n;  % MC estimate of the value of the option
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

% exact value from Black-Scholes formula
tmp=1/(sigma*sqrt(T))*log(K/s0)+(sigma/2-r/sigma)*sqrt(T);
vexact=s0*normcdf(sigma*sqrt(T)-tmp)-K*exp(-r*T)*normcdf(-tmp);

fprintf('\n theta=%8.4e, v=%8.4e, se=%8.4e, vexact=%8.4e\n',theta, v, se, vexact);

function F=fun(theta, s0, r, sigma, T, K)
F=s0*exp(-0.5*sigma^2*T+sigma*sqrt(T)*theta)*(sigma*sqrt(T)-theta)+exp(-r*T)*K*theta;



% Estimate the price of an European call option on a stock using the cross-entropy method
% The stock price is modeled by a geometric Brownian motion.

function [v, se, vexact] = option_European_CE(s0, r, sigma, T, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.


h=@(x)payoff(x,s0,r,sigma,T,K); % define a handle for the discounted payoff function

% compute the titling parameter using pilot samples
N=2000;
x=randn(N,1); tmp=h(x);
theta = sum(tmp.*x)/sum(tmp);

%theta =0; % plain MC 

y=theta+randn(n,1); % iid samples from the tilted distribution

H=h(y).*exp(-theta*y+theta^2/2); % discounted payoff x the likelihood ratio

v = sum(H)/n;  % MC estimate of the value of the option
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

% exact value from Black-Scholes formula
tmp=1/(sigma*sqrt(T))*log(K/s0)+(sigma/2-r/sigma)*sqrt(T);
vexact=s0*normcdf(sigma*sqrt(T)-tmp)-K*exp(-r*T)*normcdf(-tmp);

fprintf('\n theta=%8.4e, v=%8.4e, se=%8.4e, vexact=%8.4e\n',theta, v, se, vexact);


function h = payoff(x,s0,r, sigma, T, K)
h=max(s0*exp(-0.5*sigma^2*T+sigma*sqrt(T)*x)-exp(-r*T)*K, 0);
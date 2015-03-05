% Estimate the price of a European call option on a stock using stratified asampling
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_euro_stratification(s0, r, sigma, T, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.

k=1; % number of strata
nk=n/k; % sample size in each strata

for i=1:k
   v=rand(nk,1);  % nk samples from uniform distribution on [0,1]
   u= (i-1)/k+v/k; 
   z= norminv(u);  % nk samples from normal distribution on the i-th stratum
   x=s0 *exp(-sigma^2*T/2 +sigma*sqrt(T)*z); % discounted terminal stock price
   h=max(x-exp(-r*T)*K, 0); % discounted payoff

   vk(i) = sum(h)/nk;  % MC estimate in the i-th stratum
   sek(i) = sqrt(sum((h - vk(i)).^2)/(nk-1));  % estimate of the standard error in the i-th stratum
end

v=sum(vk)/k;
se=sqrt(sum(sek.^2/nk))/k;

fprintf('\n v=%8.4e, se=%8.4e\n',v, se);



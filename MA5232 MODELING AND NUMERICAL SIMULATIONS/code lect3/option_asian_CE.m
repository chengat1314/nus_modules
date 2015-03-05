% Estimate the price of an Asian call option on a stock using the cross-entropy method
% The stock price is modeled by a geometric Brownian motion.

function [v, se] = option_asian_CE(s0, r, sigma, T, m, K, n)

% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; m - number of monitor dates;
% K - strike price; N: number of pilot samples
% n - sample size.

t=[0: 1/m: 1]*T;  % m monitoring dates 
dt=t(2:m+1)-t(1:m);

% compute the tilting parameter using the basic cross entropy method
tmp1=zeros(m,1);
tmp2=0;
N=2000;    % number of pilot samples
for i=1:N
  z=randn(m,1);
  s=gBM(s0, r, sigma, m, dt, z);
  H=exp(-r*T)*max(sum(s(2:m+1))/m -K,0); % payoff for the k-th path

  tmp1=tmp1+H*z;
  tmp2=tmp2+H;
end
theta= tmp1/tmp2;   % tilting parameter

%theta(1:m,1)=0; %plain MC

% important sampling
for i=1:n    % loops to generate n sample paths and payoffs
  y=randn(m,1)+theta;
  s=gBM(s0,r,sigma,m,dt,y);
  H(i)=exp(-r*T)*max(sum(s(2:m+1))/m -K,0); %discounted payoff
  H(i)=H(i)*exp(-sum(theta.*y)+0.5*sum(theta.^2)); % likelihood ratio
end

v=sum(H)/n; 
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('\n v=%8.4e, se=%8.4e\n', v, se);


% generate a discrete path of geometric Bronian motion
function s=gBM(s0, r, sigma, m, dt, z)
s(1) =s0;
for k=2:m+1
s(k) = s(k-1)*exp((r-0.5*sigma^2)*dt(k-1)+sigma*sqrt(dt(k-1))*z(k-1));
end

function homework1_2()
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
ss(1) = s0;
%% Arithmetic sampling
n=n1; 
for k=1:n    % loops to generate n sample paths and payoffs
  for i=2:m+1   % generate the k-th path and payoff  
    z=randn;  % generate a sample from N(0,1)
    s(i) = s(i-1)*exp((r-0.5*sigma^2)*(t(i)-t(i-1))+sigma*sqrt(t(i)-t(i-1))*z);
    ss(i) = ss(i-1)*exp((r-0.5*sigma^2)*(t(i)-t(i-1))-sigma*sqrt(t(i)-t(i-1))*z);
  end
  x(k)=max(max(s) -K,0); % payoff for the k-th path
  y(k)=max(max(ss) -K,0);
end
v = sum(x+y)/(2*n);  % MC estimate of the value of the option
se = sqrt((sum(((x+y)/2).^2) - n*v^2)/n/(n-1));  % estimate of the standard error
fprintf('the result of Arithmetic sampling: \n v=%8.4e, se=%8.4e\n',v, se);
function homework2_2()
clear;clc;
% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.
n=10000;
a=1;
h=@(x)payoff(x,a); % define a handle for the discounted payoff function

% compute the titling parameter using pilot samples
N=2000;
x=randn(N,1); tmp=h(x)';
theta = sum(tmp.*x)/sum(tmp);

%theta =0; % plain MC 

y=theta+randn(n,1); % iid samples from the tilted distribution

H=h(y)'.*exp(-theta*y+theta^2/2); % discounted payoff x the likelihood ratio

v = sum(H)/n;  % MC estimate of the value of the option
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

fprintf('the cross-entropy method \n theta=%8.4e, v=%8.4e, se=%8.4e \n',theta, v, se);

function Hx = payoff(x,a) 
for i = 1:length(x)
   if  x(i)>0
       Hx(i) = exp(a*x(i)^(1/2));
   else 
       Hx(i) = 0 ;
   end
end 

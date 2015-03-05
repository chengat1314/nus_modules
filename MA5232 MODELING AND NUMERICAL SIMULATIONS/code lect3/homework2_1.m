function homework2_1()
clc;clear;
% mode matching to determine theta
% s0 - initial stock price; r - risk-free interest rate;
% sigma - volatility; T - terminal time; K - strike price;
% n - sample size.
a = 1;
n=10000; 
n=5*n
f=@(theta)fun(theta, a); % define a function handle f
theta = fzero(f, 0.5); % find the zero of f; or using fsolve(f, 10)

%theta =0; % plain MC 

y=theta+randn(n,1); % iid samples from the tilted distribution

H= hx(a,y)'; %discouted payoff
H=H.*exp(-theta*y+theta^2/2); % discounted payoff x the likelihood ratio

v = sum(H)/n;  % MC estimate of the value of the option
se = sqrt((sum(H.^2) - n*v^2)/n/(n-1));  % estimate of the standard error

% exact value 
% n=10*n;
randv = randn(n,1);
for i = 1:n
   if  randv(i)>0
       ex(i) = exp(a*randv(i)^(1/2));
   else 
       ex(i) = 0 ;
   end
end 
vexact=mean(ex);
see = sqrt((sum(ex.^2) - n*vexact^2)/n/(n-1)); 
fprintf('mode matching method:\n theta=%8.4e, v=%8.4e, se=%8.4e \n',theta, v, se);
fprintf('plain Monte Carlo=%8.4e, se=%8.4e \n',vexact, see);

function F=fun(theta, a)
F= exp(a*sqrt(theta)*exp(-theta^2/2))*(-theta+(a/2)*theta^(-1/2));



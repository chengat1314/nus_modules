function stochastic_demo()
clear all;close all;clc;
n = 100;a =4;b=6;xmin = 1;xmax=n; ymin =a*xmin+b;ymax=a*xmax+b;
x = [1:0.5:n];m = length(x)
y = a*x + b + (randn(1,m))*b*0.1;
plot(x,y,'.');  hold on;
xx=[ones(m,1),x']; yy=y';
w_right =  inv(xx'*xx)*xx'*yy  %best solution

p_yy=w_right(1)+w_right(2).*x; plot(x,p_yy,'-r'); hold on;
%% batch training
lambda = 500; %regularization 
[w_batch,batchthetaRecord, batchCostHistory] = gradientDescent(xx,[0;0;],yy,1e-2/n,lambda,100);
w_batch
p_yy=w_batch(1)+w_batch(2).*x; plot(x,p_yy,'-b'); hold on;
%% streaming training ;
lambda = 5000; %regularization
[w_stream,thetaRecord,  StochasticCostHistory] = StochasticGradientDescent(xx,[0;0;],yy,1e-2/n,lambda,100);
w_stream
p_yy=w_stream(1)+w_stream(2).*x; plot(x,p_yy,'-y'); hold on;
%% plot cost funciton 
figure(2)
plot([1:length(batchCostHistory)]*m,batchCostHistory);
hold on
plot(StochasticCostHistory,'r');xlabel('n');ylabel('cost')
ylim([0,100]);legend('Batch Cost','Stochastic Cost')
hold off
%% accuracy
[1-norm(w_right-w_batch)/norm(w_right),1-norm(w_right-w_stream)/norm(w_right)]
%% step in 
figure(5)
plotme(a,b,x,y)
hold on
plot(batchthetaRecord(:,2),batchthetaRecord(:,1),'*-r')
%% step in 
figure(6)
plotme(a,b,x,y)
hold on
plot(thetaRecord(:,2),thetaRecord(:,1),'*-b')

function plotme(a,b,x,y)
ma = linspace(3,6);
mb = linspace(-20,20);
[mx,my] = meshgrid(ma,mb);
for i=1:length(mx)
    for j=1:length(my)
        mz(i,j) = sum((mx(i,j).*x+my(i,j)-y).^2) ;
    end
end
contour(mx,my,mz,100) %,'ShowText','on'





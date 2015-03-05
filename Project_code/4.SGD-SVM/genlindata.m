% generate data for SVM classification/probability regression

clear

C = 10;		% cost, soft margin
N = 2;  	% dimension of X
L = 100;	% number of samples

% generate data
x = randn(L,N);
w0 = [4;-2];
b0 = 0.2;
y = 2*(tanh(x*w0+b0)>2*rand(L,1)-1)-1;	% logistic

% save data
save lindata

% plot
x0 = x(find(y<0),:);
x1 = x(find(y>0),:);
plot(x0(:,1),x0(:,2),'go',x1(:,1),x1(:,2),'g+');

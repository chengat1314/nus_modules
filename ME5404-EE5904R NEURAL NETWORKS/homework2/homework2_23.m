function homework2_23()
clc,clear,close all;
%% BP MLP batch_training_trainlm
fx = @(x)(1.2*sin(pi*x)-cos(2.4*pi*x));
x = [-1:0.01:1];
% plot(x,fx(x),'r.')
%% Generate the training data
xtrain = [-1:0.05:1];
xtest = [-1:0.01:1];
P = xtrain; T = fx(xtrain);  
PP = xtest; TT = fx(xtest); 
sizen=[1:10,50];
mse = zeros(2,length(sizen));
testx = zeros(2,length(sizen));
%% BP batch training trainbr test
for n = 1:length(sizen)
    figure(n)
    [net,net_initial] = batch_training_trainbr(P,T,n);
    initial_Y = sim(net_initial,P);  % 
    plot(P,T,P,initial_Y,'.') ; hold on
    Y = sim(net,P);  % training error  
    plot(P,Y,'*r')
    hold on
    YY = sim(net,PP);  % testing error  
    plot(PP,YY,'.-b')
    title('batch_training_trainbr');
    xlabel('x');ylabel('f(x)');
    legend('the orginal line','initial point','training result','test samples');
    train_error = sum((T-Y).^2)/(length(Y)-2);
    test_error = sum((TT-YY).^2)/(length(YY)-2);
    mse(:,n) = [train_error;test_error];
    YX = sim(net,[-1.5,1.5]); 
    testx(:,n) = YX;
end
figure(20)
plot(sizen,mse(1,:),'*-');hold on
plot(sizen,mse(2,:),'o-');legend('train error','test error');
xlabel('n');ylabel('mse');
[testx',testx(1,:)' - fx(-1.5),testx(2,:)' - fx(1.5)]




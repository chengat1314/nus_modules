function [net,net_initial] = batch_training_trainbr(P,T,n)
net = newff(P,T,5,{},'trainbr'); %% create s feed back progation network
net_initial = net;
net.trainParam.epochs = 1000; 
net = train(net,P,T);  
%% Description
% trainbr is a network training function that updates the weight and bias values 
% according to Levenberg-Marquardt optimization. 
% It minimizes a combination of squared errors and weights, 
% and then determines the correct combination so as to produce a network 
% that generalizes well. The process is called Bayesian regularization.
function [net,net_initial] = batch_training_trainlm(P,T,n)
net = newff(P,T,5,{},'trainlm'); %% create s feed back progation network
net_initial = net;
net.trainParam.epochs = 1000; 
net = train(net,P,T);

%% Description
% trainlm is a network training function that updates weight 
% and bias values according to Levenberg-Marquardt optimization.
% 
% trainlm is often the fastest backpropagation algorithm in the toolbox,
% and is highly recommended as a first-choice supervised algorithm, 
% although it does require more memory than other algorithms.
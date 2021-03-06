function [net,net_initial] = sequential_training(P,T,n)
%% P input; T output; n number of the neuron
net = newff(P,T,n,{'tansig' 'purelin'},'trainlm');
net_initial = net;
epoch_s = 1000;  % Specify the number of epoch for sequential training
net.trainParam.epochs = 1;  % Set epoch of batch training to be 1

for i = 1 : epoch_s
    index = randperm(length(P));   % Shuffle the input data every epoch
    net = adapt(net,P(index),T(index));    % Perform sequential learning
end

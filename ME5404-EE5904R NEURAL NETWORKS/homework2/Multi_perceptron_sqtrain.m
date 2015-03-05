function error = Multi_perceptron_sqtrain(n)
[trainB1,trainB2,testB1,testB2]= readimg();
% size(trainB1),size(trainB2),size(testB1),size(testB2)

P = [trainB1,trainB2];
T = [ones(1,50),zeros(1,50)];
PP = [testB1,testB2];
TT = [ones(1,15),zeros(1,15)];
error = zeros(2,n);

epoch_s = 100;  % Specify the number of epoch for sequential training
net.trainParam.epochs = 1;  % Set epoch of batch training to be 1
for i = 1:n
    net = newff(P,T,i,{'tansig' 'purelin'},'trainlm');
    %     view(net);close all
    for j = 1 : epoch_s
        index = randperm(100); % Shuffle the input data every epoch
        net = adapt(net,P(:,index),T(:,index));    % Perform sequential learning
    end
    trainresult = sim(net,P);
    trainer = threshold(trainresult);
    trainerror = 1 - sum(trainer==T)/length(trainer);
    %% test
    
    testresult = sim(net,PP);
    tester = threshold(testresult);
    testerror = 1 - sum(tester==TT)/length(tester);
    error(:,i) = [trainerror;testerror];
end
function error = single_perceptron()
%% single perceptron
[trainB1,trainB2,testB1,testB2]= readimg();
% size(trainB1),size(trainB2),size(testB1),size(testB2)

P = [trainB1,trainB2];
T = [ones(1,50),zeros(1,50)];
net = newff(P,T,1,{'purelin'},'traincgf');
net.trainParam.epochs = 1000; 
net = train(net,P,T);
% view(net);close all
trainresult = sim(net,P);
trainer = threshold(trainresult);
trainerr = 1- sum(trainer==T)/length(trainer);
%% test 
PP = [testB1,testB2];
TT = [ones(1,15),zeros(1,15)];
testresult = sim(net,PP);
tester = threshold(testresult);
testerr = 1- sum(tester==TT)/length(tester);
error =[trainerr,testerr];
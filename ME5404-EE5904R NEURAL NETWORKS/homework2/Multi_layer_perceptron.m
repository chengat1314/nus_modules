function error = Multi_layer_perceptron(n)

[trainB1,trainB2,testB1,testB2]= readimg();
% size(trainB1),size(trainB2),size(testB1),size(testB2)

P = [trainB1,trainB2];
T = [ones(1,50),zeros(1,50)];
PP = [testB1,testB2];
TT = [ones(1,15),zeros(1,15)];
error = zeros(2,n);
for i = 1:n
    net = newff(P,T,i,{'tansig' 'purelin'},'traincgf');
    net.trainParam.epochs = 1000; 
    net = train(net,P,T);
%     view(net);close all
    trainresult = sim(net,P);
    trainer = threshold(trainresult);
    trainerror = 1 - sum(trainer==T)/length(trainer);
    %% test 

    testresult = sim(net,PP);
    tester = threshold(testresult);
    testerror = 1 - sum(tester==TT)/length(tester);
    error(:,i) = [trainerror;testerror];
end
% plot(error(1,:),'*-');hold on
% plot(error(2,:),'o-');legend('train error','test error');
% xlabel('n');ylabel('mse');title('training error plot');
% error
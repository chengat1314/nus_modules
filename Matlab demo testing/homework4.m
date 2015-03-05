function homework4()
clc;clear;
train_data=load('D:/hw4_nnet_train.dat');
test_data=load('D:/hw4_nnet_test.dat');
[m1,n1]=size(train_data);
[m2,n2]=size(test_data);
%% plot
% data1 = train_data(train_data(:,3)==1,:);
% data2 = train_data(train_data(:,3)==-1,:);
% plot(data1(:,1),data1(:,2),'*r');
% hold on
% plot(data2(:,1),data2(:,2),'ob');
% data1 = test_data(test_data(:,3)==1,:);
% data2 = test_data(test_data(:,3)==-1,:);
% figure(2)
% plot(data1(:,1),data1(:,2),'*r');
% hold on
% plot(data2(:,1),data2(:,2),'ob');
%%
train_data(train_data(:,3)==-1,3)=0;
test_data(test_data(:,3)==-1,3)=0;
P = train_data(:,1:2)'; T = train_data(:,3)';
PP = test_data(:,1:2)'; TT = test_data(:,3)';
m = 1;
ni = [1,6,11,16,21];
ri=[0.001,0.01,0.1,1,10];
record = zeros(m,length(ri));
%% create s feed back progation network trainbr with regularization
for k = 1:m
    e_out = zeros(1,length(ri));
    for i = 1:length(ri)
        %net = newff(P,T,3,{'tansig'},'trainbr'); 
        net = newff(P,T,[8 3 1],{'tansig','tansig','tansig','tansig'},'trainbr'); 
%        net.trainParam.lr = 0.01; %learning_rate
%         net.trainParam.mu_max =  0.1;
 %       net.performParam.regularization = 0.1;
        % uniformly from the range (-0.1,0.1)
        net.trainParam.epochs = 1000; 
%         net.trainParam.time = 100;
        net = train(net,P,T);
        YY = threshold_01(sim(net,PP));   % testing error  
        e_out(i) = 1 - sum(YY==TT)*1.0/length(YY);
    end
    record(k,:)=e_out;
    fprintf('times: %i\n',k)
end
mean(record)


%0.6631    0.6494    0.5889    0.6076    0.6271
%0.6588    0.6540    0.6271    0.5660    0.6205
% 0.6861    0.6690    0.6178    0.5861
% 0.6346    0.6546    0.6939    0.6215    0.7049


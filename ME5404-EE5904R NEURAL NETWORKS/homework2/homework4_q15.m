function homework4_q15()
clc,clear;
test = load('D:/hw4_knn_test.dat');
train = load('D:/hw4_knn_train.dat');
y= train(:,end);
yy=test(:,end);
[m2,n2] = size(test);
[m1,n1] = size(train);
%% distance
dist1 = distij(train(:,1:end-1),train(:,1:end-1));
dist2 = distij(test(:,1:end-1),train(:,1:end-1));

%% predict
% output1 = zeros(m1,1);
% for i=1:m1
%     out = mindist(dist1(i,:),i);
%     output1(i) = y(out);
% end
% ein = 1-(sum(output1==y)/length(output1))
% output2 = zeros(m2,1);
% for i=1:m2
%     [v,out] = min(dist2(i,:));
%     output2(i) = y(out);
% end
% eout = 1-(sum(output2==yy)/length(output2))

%% knn predict
output1 = knnclassify(train(:,1:end-1), train(:,1:end-1), y, 5); 

ein = 1-(sum(output1==y)/length(output1))
output2 =  knnclassify(test(:,1:end-1), train(:,1:end-1), y, 5); 
eout = 1-(sum(output2==yy)/length(output2))


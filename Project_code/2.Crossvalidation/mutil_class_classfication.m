function mutil_class_classfication()
clc,clear,close all;
load fisheriris                              %# load iris dataset
p = 0.5; group = grp2idx(species) ;
[train,test] = crossvalind('Holdout',group,p); %split the data set into two parts
trainingSample = meas(train,3:4) ;
trainingLabel = group(train,1);
testingSample = meas(test,3:4) ;
testingLabel = group(test,1);

numClass = max(trainingLabel);
DecisionMatrix = zeros(numClass,sum(test));
%multiclass training 1 visa 1
for i = 1:numClass
    for j = i+1:numClass
        indexij = (trainingLabel==i)|(trainingLabel==j);
        traingsampleij = trainingSample(indexij,:);
        trainglabelij = trainingLabel(indexij,:);
        
        k = 5; %k-fold CV
        cvFolds = crossvalind('Kfold', trainglabelij, k);  %# get indices of k-fold CV
        sigma = 2.^[-5:5];
        C = 2.^[-5:5];
        [bestsigma,bestC] = selectBestParametersRBF(traingsampleij,trainglabelij,sigma,C,cvFolds,k)
%         figure(i+j)
        svmModel = svmtrain(traingsampleij, trainglabelij, ...
            'Showplot',0, 'Method','QP', ...
            'BoxConstraint',bestC,'Kernel_Function','rbf', 'RBF_Sigma',bestsigma);
        
%         title(['Para-Selection by Cross-validation the best C=' num2str(bestC) ', sigma=' num2str(bestsigma)]);
        
        pred = svmclassify(svmModel, testingSample, 'Showplot',false);
        %acc(i,j) = sum(grp2idx(pred)==grp2idx(trainingLabel))/length(trainingLabel)
        rowIndex = pred;
        colIndex = [1:sum(test)]';
        linearInd = sub2ind([numClass sum(test)],rowIndex,colIndex);
        DecisionMatrix(linearInd) = DecisionMatrix(linearInd) + 1;
    end
end

[temp,outlabel] = max(DecisionMatrix);
acc = mean(grp2idx(outlabel) == grp2idx(testingLabel)) 
%three classes prediction accuracy



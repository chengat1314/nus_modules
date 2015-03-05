function  svm_crossvalidation()
clc,clear,close all;
load fisheriris                              %# load iris dataset 
trainingSample = meas(51:end,3:4);
trainingLabel = species(51:end,1);
k = 5; %k-fold CV 
cvFolds = crossvalind('Kfold', trainingLabel, k);  %# get indices of k-fold CV
sigma = 2.^[-5:5]; 
C = 2.^[-5:5]; 
[bestsigma,bestC]=selectBestParametersRBF(trainingSample,trainingLabel,sigma,C,cvFolds,k)

svmModel = svmtrain(trainingSample, trainingLabel, ...
    'Showplot',1, 'Method','QP', ... 
    'BoxConstraint',bestC,'Kernel_Function','rbf', 'RBF_Sigma',bestsigma);
title(['Para-Selection by Cross-validation the best C=' num2str(bestC) ', sigma=' num2str(bestsigma)]);

pred = svmclassify(svmModel, trainingSample, 'Showplot',false);
acc = sum(grp2idx(pred)==grp2idx(trainingLabel))/length(trainingLabel)


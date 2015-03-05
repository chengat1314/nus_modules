function [bestsigma,bestC]=selectBestParametersRBF(trainingSample,trainingLabel,sigma,C,cvFolds,k) 
accurcy = zeros(length(C),length(sigma));
for sc = 1:length(C)
    for  j = 1:length(sigma)
%         outlabel = zeros(length(trainingLabel),1);
        for i = 1:k                                 %# for each fold
            testingFoldSample = trainingSample(cvFolds == i,:);  %# get indices of test instances
            trainingFoldSample =  trainingSample(cvFolds ~= i,:);                    %# get indices training instances
            trainingFoldLabel =  trainingLabel(cvFolds ~= i,:);                    %# get indices training instances

            %# train an SVM model over training instances
            svmModel = svmtrain(trainingFoldSample, trainingFoldLabel, ...
                          'Showplot',false, 'Method','QP', ...  %'Autoscale',true,'BoxConstraint',2e-1,
                          'BoxConstraint',C(sc),'Kernel_Function','rbf', 'RBF_Sigma',sigma(j));

            %# test using test instances
            pred = svmclassify(svmModel, testingFoldSample, 'Showplot',false);
            outlabel(cvFolds == i) = pred;
            %# evaluate and update performance object
            % cp = classperf(cp, pred, testIdx);
        end
        %# get accuracy
        accurcy(sc,j) = sum(grp2idx(outlabel)==grp2idx(trainingLabel))/length(outlabel); %cp.CorrectRate;
        %# get confusion matrix  columns:actual, rows:predicted, last-row: unclassified instances
    %     cp.CountingMatrix
    end
end
% accurcy ;%two dimensional accurcy table
[maxCol,Icol]= max(accurcy);%find the maximun accurcy
[acc,IRow] = max(maxCol);
bestsigma = sigma(IRow); %
bestC = C(Icol(IRow));


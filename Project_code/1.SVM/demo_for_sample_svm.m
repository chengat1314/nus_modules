function demo_for_sample_svm()
clc,clear;
load fisheriris                              %# load iris dataset
% groups = ismember(species,'setosa');         %# create a two-class problem
v=[4,3];
x = meas(:,v);%1-4, choose variable
gscatter(x(:,1),x(:,2),species);

%% firstly consider 2d categories
close all
x = meas(51:end,v);
y = species(51:end,1);
gscatter(x(:,1),x(:,2),y);

%% linear SVM
figure(2)
subplot(1,3,1)
linearsvm = svmtrain(x,y,'showplot',true);
title('linear kernel')
%% nonlinear polynomial SVM
% figure(3)
subplot(1,3,2)
nonlinear = svmtrain(x,y,'showplot',true,'kernel_function','polynomial');
title('polynomial kernel n=3')
% figure(6)
subplot(1,3,3)
nonlinear = svmtrain(x,y,'showplot',true,'kernel_function','polynomial','Polyorder',10);
title('polynomial kernel n=10')
%% nonlinear RBF SVM
figure(4)
subplot(1,3,1)
nonlinear = svmtrain(x,y,'showplot',true,'kernel_function','rbf','rbf_sigma',0.1);
title('RBF kernel sigma=0.1')
subplot(1,3,2)
nonlinear = svmtrain(x,y,'showplot',true,'kernel_function','rbf','rbf_sigma',1);
title('RBF kernel sigma=1')
subplot(1,3,3)
nonlinear = svmtrain(x,y,'showplot',true,'kernel_function','rbf','rbf_sigma',10);
title('RBF kernel sigma=10')



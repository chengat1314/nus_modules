function [ypred,indw] = svcm_test(xtest,ytest,xtrain,ytrain,atrain,btrain);
% function [ypred,indw] = svcm_test(xtest,ytest,xtrain,ytrain,atrain,btrain);
%
%	support vector classification machine
%	test set cross-validation
%	soft margin
%	uses "kernel.m"
%
%	xtrain: (Ltrain,N) with Ltrain: number of points; N: dimension
%	ytrain: (Ltrain,1) containing class labels (-1 or +1)
%	xtest:  (Ltest,N) with Ltest: number of points; N: dimension
%	ytest:  (Ltest,1) containing class labels (-1 or +1)
%	atrain: alpha coefficients (from svcm_train on xtrain and ytrain)
%	btrain: offest coefficient (from svcm_train on xtrain and ytrain)
%
%	ypred:  predicted y; (Ltest,1) containing class labels (-1 or +1)
%	indw:   indices of wrongly classified test points

[Lxtrain,Nxtrain] = size(xtrain);
[Lytrain,Nytrain] = size(ytrain);
[Lxtest,Nxtest] = size(xtest);
[Lytest,Nytest] = size(ytest);
if Nxtrain~=Nxtest
    fprintf('svcm_test error: xtrain and xtest different number of dimensions (%g/%g)\n',...
	 Nxtrain, Nxtest);
    return
elseif Lytrain~=Lxtrain
    fprintf('svcm_test error: xtrain and ytrain different number of data points (%g/%g)\n',...
	 Lxtrain, Lytrain);
    return
elseif Lytest~=Lxtest
    fprintf('svcm_test error: xtest and ytest different number of data points (%g/%g)\n',...
	 Lxtest, Lytest);
    return
elseif Nytrain~=1
    fprintf('svcm_test error: ytrain not a single variable (%g)\n', Nytrain);
    return
elseif Nytest~=1
    fprintf('svcm_test error: ytest not a single variable (%g)\n', Nytest);
    return
elseif any(ytrain~=-1&ytrain~=1)
    fprintf('svcm_test error: ytrain takes values different from {-1,+1}\n');
    return
elseif any(ytest~=-1&ytest~=1)
    fprintf('svcm_test error: ytest takes values different from {-1,+1}\n');
    return
end

fprintf('Support vector soft-margin classifier\n')
fprintf('%4g test points\n', Lxtest)
fprintf('%4g dimensions\n', Nxtest)

ac = atrain(atrain>0);
xc = xtrain(atrain>0,:);
yc = ytrain(atrain>0);

ypred = 2*(kernel(xtest,xc)*(ac.*yc)+btrain>0)-1;	% +/- 1
indw = find(ypred~=ytest);

nw = length(indw);
fprintf('%4g test points incorrectly classified (%3.1f%%)\n\n', nw, nw/Lxtest*100)

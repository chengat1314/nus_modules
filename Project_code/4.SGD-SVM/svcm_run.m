function [ypred,margin] = svcm_run(xrun,xtrain,ytrain,atrain,btrain);
% function [ypred,margin] = svcm_run(xrun,xtrain,ytrain,atrain,btrain);
%
%	support vector classification machine
%	soft margin
%	uses "kernel.m"
%
%	xtrain: (Ltrain,N) with Ltrain: number of points; N: dimension
%	ytrain: (Ltrain,1) containing class labels (-1 or +1)
%	xrun:   (Lrun,N) with Lrun: number of points; N: dimension
%	atrain: alpha coefficients (from svcm_train on xtrain and ytrain)
%	btrain: offest coefficient (from svcm_train on xtrain and ytrain)
%
%	ypred:  predicted y; (Lrun,1) containing class labels (-1 or +1)
%	margin: (signed) separation from the separating hyperplane; (Lrun,1)

[Lxtrain,Nxtrain] = size(xtrain);
[Lytrain,Nytrain] = size(ytrain);
[Lxrun,Nxrun] = size(xrun);
if Nxtrain~=Nxrun
    fprintf('svcm_run error: xtrain and xrun different number of dimensions (%g/%g)\n',...
	 Nxtrain, Nxrun);
    return
elseif Lytrain~=Lxtrain
    fprintf('svcm_run error: xtrain and ytrain different number of data points (%g/%g)\n',...
	 Lxtrain, Lytrain);
    return
elseif Nytrain~=1
    fprintf('svcm_run error: ytrain not a single variable (%g)\n', Nytrain);
    return
elseif any(ytrain~=-1&ytrain~=1)
    fprintf('svcm_run error: ytrain takes values different from {-1,+1}\n');
    return
end

fprintf('Support vector soft-margin classifier\n')
fprintf('  %g run points\n', Lxrun)
fprintf('  %g dimensions\n\n', Nxrun)

nonzero = find(atrain~=0);
ac = atrain(nonzero);
xc = xtrain(nonzero,:);
yc = ytrain(nonzero);

margin = kernel(xrun,xc)*(ac.*yc)+btrain;   % separation from the separating hyperplane
ypred = 2*(margin>0)-1;                     % thresholded, +/- 1

fprintf('  (%g/%g) run points classified as (-1/+1)\n\n', sum(ypred==-1), sum(ypred==1))

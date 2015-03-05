clear all

global KTYPE
global KSCALE

KTYPE = 6;	% gaussian kernel
KSCALE = .025;	% 1/(2*sigma)^2 for gaussian kernel
C = 1;		% soft margin regularization parameter

load diabetes
[L,N] = size(x);
if L~=768 | N~=8
    fprintf('diabetes data error: wrong matrix dimensions\n')
    return
end
Ltrain = 576;

xtrain = x_norm(1:Ltrain,:);
ytrain = y(1:Ltrain);
xtest =x_norm(Ltrain+1:L,:);
ytest =y(Ltrain+1:L);

[a,b,D,inds,inde,indwLOO] = svcm_train(xtrain,ytrain,C);
[ypred,indwTEST] = svcm_test(xtest,ytest,xtrain,ytrain,a,b);

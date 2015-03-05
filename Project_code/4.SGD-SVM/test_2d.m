clear all

global KTYPE
global KSCALE
global online
global visualize

KTYPE = 6;	% gaussian kernel
KSCALE = .25;	% 1/(2*sigma)^2 for gaussian kernel
online = 1
visualize = 1;

C = 10;		% regularization parameter
N = 2;          % dimension of X; 2!
L = 100;        % number of datapoints


% generate or load 2D data

% linear data
% w0 = [4;-2];
% b0 = 0.2;
% [x,y] = genlindata(L,N,w0,b0);          % training set

% nonlinear data
% L0 = round(L/4)*ones(4,1);
% x0 = [0.9,1.5;-1.1,0.8;-0.7,-1.3;1.4,-0.5];
% s0 = [1.3;0.8;1.1;0.9];
% y0 = [1;-1;1;-1];
% [x,y] = gennonlindata(N,L0,x0,s0,y0);     % training set

% backup data
% save data

% retrieve data
load data

% train SVM

[a,b,D,inds,inde,indwLOO] = svcm_train(x,y,C);

figure(1)
SetFont('Palatino', 18, 16, 18, 15);
SetMarker(2);
print -deps gctraj.eps

figure(2)
SetFont('Palatino', 18, 16, 18, 15);
print -deps atraj.eps

figure(3)
SetFont('Palatino', 18, 16, 18, 15);
print -depsc gtraj.eps  % color


% run SVM on test grid

xmax = max(max(abs(x)));
nx = 50;
ngrid = 2*nx+1; % 101 by 101 mesh
ngrid2 = ngrid^2;
xgrid = (-nx:nx)'*xmax/nx;
grid = ones(ngrid,1)*xgrid';
xtest(:,1) = reshape(grid,ngrid2,1);
xtest(:,2) = reshape(grid',ngrid2,1);
[ypred,margin] = svcm_run(xtest,x,y,a,b);
mgrid = reshape(margin,ngrid,ngrid);
mmax = max(abs(margin));

ind0=find(y==-1);
x0=x(ind0,:);                 % zeros
ind1=find(y==1);
x1=x(ind1,:);                 % ones
inds0=intersect(inds,ind0);
xs0=x(inds0,:);               % support vectors, zeros
inds1=intersect(inds,ind1);
xs1=x(inds1,:);               % support vectors, ones
inde0=intersect(inde,ind0);
xe0=x(inde0,:);               % error vectors, zeros
inde1=intersect(inde,ind1);
xe1=x(inde1,:);               % error vectors, ones

figure(4)
clf
grey = gray;                            
redblue = min(grey,1-grey);             
redblue(:,1) = redblue(:,1)+grey(:,1);  
redblue(:,2) = 2*redblue(:,2);          
redblue(:,3) = redblue(:,3)+1-grey(:,3);
colormap(redblue)

image(xgrid,xgrid,(0.15*mgrid+mmax)/mmax*0.5*length(gray))
colormap(redblue)
hold on
plot(x0(:,1),x0(:,2),'bv')
plot(x1(:,1),x1(:,2),'r^')
SetMarker(3)
plot(xs0(:,1),xs0(:,2),'bo')
plot(xs1(:,1),xs1(:,2),'ro')
plot(xe0(:,1),xe0(:,2),'bx')
plot(xe1(:,1),xe1(:,2),'rx')
contour(xgrid,xgrid,mgrid,[-1 -1],'b--')
contour(xgrid,xgrid,mgrid,[0 0],'k-')
contour(xgrid,xgrid,mgrid,[1 1],'r--')
hold off

axis(1.05*[-xmax,xmax,-xmax,xmax])
axis square
xlabel('{\it{x}}_1')
ylabel('{\it{x}}_2')
SetFont('Palatino', 18, 16, 18, 15);
set(gca,'YTick',get(gca,'XTick'))

for i=1:L
    h = text(x(i,1), x(i,2), [' ',int2str(i)]);
    set(h, 'FontName', 'Helvetica')
    set(h, 'FontSize', 8)
end

print -deps points.eps
print -depsc pointsc.eps % color

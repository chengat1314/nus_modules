function outlier_detection()

%% Exclude Data Using a Simple Rule
% [x, y] = titanium;
% f1 = fit(x',y','gauss2', 'Exclude', x<800);
% plot(f1,x,y,x<800)

%%  Exclude Data by Distance from the Model
xdata = (0:0.1:2*pi)';
y0 = sin(xdata);
% Response-dependent Gaussian noise
gnoise = y0.*randn(size(y0));

% Salt-and-pepper noise
spnoise = zeros(size(y0));
p = randperm(length(y0));
sppoints = p(1:round(length(p)/5));
spnoise(sppoints) = 5*sign(y0(sppoints));

ydata = y0 + gnoise + spnoise;
f = fittype('a*sin(b*x)');
fit1 = fit(xdata,ydata,f,'StartPoint',[1 1]);

fdata = feval(fit1,xdata);
I = abs(fdata - ydata) > 1.5*std(ydata);
outliers = excludedata(xdata,ydata,'indices',I);

fit2 = fit(xdata,ydata,f,'StartPoint',[1 1],...
           'Exclude',outliers);
fit3 = fit(xdata,ydata,f,'StartPoint',[1 1],'Robust','on');

plot(fit1,'r-',xdata,ydata,'k.',outliers,'m*')
hold on
plot(fit2,'c--')
plot(fit3,'b:')
xlim([0 2*pi])


       
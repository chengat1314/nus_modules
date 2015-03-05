function reg_plot()
clear;clc
% file = 'E:/courses/Course 2014 semester 2/ST5202 Applied Regression Analysis/mid-exam/cars.txt'
% fid   = fopen('E:/courses/Course 2014 semester 2/ST5202 Applied Regression Analysis/mid-exam/cars.txt');
% % data = dlmread('');
% line1 = textscan(fid, '%s%s%s%f%f%f%f%f\r\n %*[^\n]','HeaderLines',0);
%% "CityMPG"	"Weight"	"Horsepower"
data= load('E:\courses\Course 2014 semester 2\ST5202 Applied Regression Analysis\mid-exam\data.txt');
y = 100./data(:,1);
x1 = data(:,2);
x2 = data(:,3)./data(:,2);
XX = [ones(length(x1),1),x1,x2];
 [b,bint] = regress(y,XX)
% lm = fitlm(XX,y,'linear')
% anova(lm,'summary')

plot3(x1,x2,y,'ro');
grid on
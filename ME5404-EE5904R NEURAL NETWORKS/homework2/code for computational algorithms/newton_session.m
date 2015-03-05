% Sample MATLAB session with Newton's method.
% Download codes
% newton.m and myfunc.m into your directory
% from where your run MATLAB.
%
clc
clear all;
[x,xiter,ithist,iflag] = newton(@myfunc,[0 3]',1e-6,100);
'the solution x is'
x
'the iteration history is'
ithist
'iflag is'
iflag
'x value in the 7th iteration is'
xiter{7}
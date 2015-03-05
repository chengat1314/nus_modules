function demo()
%% athour:chengf@nus
% demo for test fft2d and ifft2d
close all; clear all;clc
%% import date  source
img   = imread('D:/data/satellite.png','png');
img = double(img); %enable the matrix can caculate the complex number.
[m,n] = size(img);
M=ceil(log2(m)); N=ceil(log2(n));  %nextpow2(m)                 
%% initial the matrix
fft_x=zeros(m,n);
ifft_x=zeros(m,n);
%% test the my function 
fprintf('FFT time cost for my function:')
tic
fft_x=fft2d(img,M,N);
toc
fprintf('IFFT time cost for my function:')
tic
ifft_x=ifft2d(fft_x,M,N);
toc
%% compare with the system function 
fprintf('FFT time cost for system function:')
tic
y=fft2(img);
toc
fprintf('IFFT time cost for system function:')
tic
yy=ifft2(y);
toc
%% show the result of my fft
figure(1);
imagesc(100*log(1+abs(fftshift(fft_x)))); colormap(gray); 
title('magnitude spectrum');

figure(2);
imagesc(angle(fft_x)); colormap(gray);
title('phase spectrum');

figure(3);
imagesc(abs(ifft_x)); colormap(gray);
title('ifft image');
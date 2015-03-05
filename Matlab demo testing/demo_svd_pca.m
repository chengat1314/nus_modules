function demo_svd_pca()
%% athour:chengf@nus
% demo for test fft2d and ifft2d
close all; clear all;clc
%% import date  source
img   = imread('D:/data/satellite.png','png');
img = double(img); %enable the matrix can caculate the complex number.
[m,n] = size(img)
imagesc(img); colormap(gray); 
[s v d]=svd(img);
k = 0;
for i = 1:m
    if v(i,i)<3000
        v(i,i)=0;
        k = k+1;
    end
end
k
ss = [s(:,1:n-k),zeros(m,k)];
vv = v;
dd = [d(:,1:n-k),zeros(m,k)];
img_new = ss*vv*dd';
figure(2)
imagesc(img_new); colormap(gray); 

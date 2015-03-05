function demo()
%% find the translation between the registration image I,T

% D:/data/image_of_template.png
% D:/data/image_for_alignment.png
% image_of_template.png: Template image 
% image_for_alignment.png: Image for registration
clc,clear
% http://en.wikipedia.org/wiki/Phase_correlation
img1   = double(imread('D:/data/image_of_template.png','png'));
img2   = double(imread('D:/data/image_for_alignment.png','png'));%satellite

shift = imshift(img1,img2)

figure(2)
imagesc(img1); colormap(gray);
figure(3)
imagesc(img2); colormap(gray); 

function shift = imshift(I,T)
%% find the translation between the registration image I,T
% img3   = imread('D:/data/satellite.png','png'); 
I = double(I);
T = double(T);
fft2d_I = fft2(I);
fft2d_T = fft2(T);

%% define the cross power spectrum
R = ((fft2d_I.*conj(fft2d_T))./abs((fft2d_I.*conj(fft2d_T))));
ifft2d_R = ifft2(double(R)) ;
abs_ifft2d_R = abs(ifft2d_R);

%% plot the delta(x-u,y-v) function
imagesc(abs_ifft2d_R); colormap(gray);

%% find the maximum point and caculate the location
[value, location] =  max(abs_ifft2d_R(:));
[m,n]=size(abs_ifft2d_R);
u = floor(location/m);
v = location - u*m;
shift = [u,v];
 
end

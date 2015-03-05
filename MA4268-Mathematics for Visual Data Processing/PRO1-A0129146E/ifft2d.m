function x = ifft2d(img,M,N)
%% author:chengf @nus 20141006
%inverse fast foureir transform with 2D data, like as image processing
%% data sources
% close all; clear all;clc
% img   = imread('D:/data/satellite.png','png');
% img = double(img);
% img = fft2(img);
[m,n] = size(img);
% M=ceil(log2(m));
% N=ceil(log2(n));
%% mark up the missed dimension 2^i columns
% tic
x1=zeros(2^M-m,n);
x=cat(1,img,x1);   % cat(1,A,B) is the same as [A;B].
x2=zeros(2^M,2^N-n);
x=cat(2,img,x2);   %  cat(2,A,B) is the same as [A,B].
%% change the position of the column
mp1=0:2^M-1;
mb=zeros(1,2^M);
for t=1:M
    mp2=floor(mp1/2);
    mb=mb*2+(mp1-2*mp2);
    mp1=mp2;
end
x(:,:)=x(mb+1,:);

np1=0:2^N-1;
nb=zeros(1,2^N);
for t=1:N
    np2=floor(np1/2);
    nb=nb*2+(np1-2*np2);
    np1=np2;
end
x(:,:)=x(:,nb+1);
%% caculate the exp iterms first
Max = max(M,N);
t=0:2^(Max-1)-1;
NN = 2^Max;
iexp(t+1)=exp(2*pi*i*t/NN); %1/N ()

% tic
for r=1:M
    B=2^(r-1); %
    for s=0:B-1
        P=s*2^(M-r);
        for k=s+1:2^r:2^M
            temp=iexp(P+1)*(x(k+B,:));
            x(k+B,:)=x(k,:) - temp;
            x(k,:)=x(k,:) + temp;
        end
    end
end

for r=1:N
    B=2^(r-1);
    for s=0:B-1
        P=s*2^(N-r);
        for k=s+1:2^r:2^N
            temp = iexp(P+1) * (x(:,k+B));
            x(:,k+B) = x(:,k) - temp;
            x(:,k) = x(:,k) + temp;
        end
    end
end
% toc
x = (1.0/NN)*x;


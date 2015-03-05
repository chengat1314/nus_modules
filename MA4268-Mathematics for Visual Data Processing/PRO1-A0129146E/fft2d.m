function x=fft2d(img,M,N)
%% author:chengf @nus 20141006
%fast foureir transform with 2D data, like as image processing
%% data sources
% close all; clear all;clc
% img   = imread('C:/Users/A0129146/Documents/course/Projects/code/satellite.png','png');
% %  img = [1:8;2:9;3:10;4:11;5:12;6:13;7:14;8:15];
% img = double(img);
% tic
[m,n] = size(img);
% M=ceil(log2(m));                  
% N=ceil(log2(n)); 
%% mark up the missed dimension 2^i columns
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
%mb = 0     4     2     6     1     5     3     7
x(:,:)=x(mb+1,:);

np1=0:2^N-1;   
nb=zeros(1,2^N);
for t=1:N
    np2=floor(np1/2);
    nb=nb*2+(np1-2*np2);
    np1=np2;
end
x(:,:)=x(:,nb+1); 

Max = max(M,N);
t=0:2^(Max-1)-1;                     
iexp(t+1)=exp(-2*pi*i*t/2^Max); 

for r=1:Max                         
    B=2^(r-1); %
    for s=0:B-1
        P=s*2^(Max-r);
        for k=s+1:2^r:2^Max   
            temp1=iexp(P+1)*(x(k+B,:));
            x(k+B,:)=x(k,:) - temp1;
            x(k,:)=x(k,:) + temp1;
            temp2 = iexp(P+1) * (x(:,k+B));
            x(:,k+B) = x(:,k) - temp2;
            x(:,k) = x(:,k) + temp2;
        end
    end 
end   


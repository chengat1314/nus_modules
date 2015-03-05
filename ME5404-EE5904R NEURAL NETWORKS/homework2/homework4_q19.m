function homework4_q19()
%% kmeans
clc,clear;
data = load('D:/hw4_kmeans_train.dat');
[m,n]=size(data);
k = 10;
[idx,C] = kmeans(data,k);
mean(data(idx==1,:));

mean(data(idx==2,:));
error = zeros(1,k);
for i=1:k
    error(i) = dist_point_center(data(idx==1,:),mean(data(idx==1,:)));
end

error
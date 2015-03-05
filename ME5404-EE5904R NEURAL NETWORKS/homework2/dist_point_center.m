function error  = dist_point_center(data ,mean)
[m,n]=size(data);
error = 0;
for i = 1:m
    error = error + sum((data(i,:)-mean).^2);
end
error =error/m;
function out = maxlike(dist1,i,k);
n = length(dist1);
temp = dist1(1);
k = 1;
for j = 1:dist1
    if(j~=i)
        if(dist1(j)<temp)
            temp = dist1(j);
            k = j;
        end
    end
end
out = k;
function output = geomean(x)
n = length(x);
a=1;
for i=1:n
    a=a*x(i);
end
output = a^(1/n);
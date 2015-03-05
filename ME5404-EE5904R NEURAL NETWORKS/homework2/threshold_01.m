function y = threshold_01(x)
for i=1:length(x)
    if x(i)>0
        y(i) = 1;
    else
        y(i) = -1;
    end
end

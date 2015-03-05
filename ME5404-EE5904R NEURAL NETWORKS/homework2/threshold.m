function y = threshold(x)
for i=1:length(x)
    if x(i)>0.5
        y(i) = 1;
    else
        y(i) = 0;
    end
end


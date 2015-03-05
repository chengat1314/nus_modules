function Hx = hx(a,x)
for i = 1:length(x)
   if  x(i)>0
       Hx(i) = exp(a*x(i)^(1/2));
   else 
       Hx(i) = 0 ;
   end
end 
function convergence_sqrt()
a=0.111;%115;%
intial = [0.01:0.01:0.6];
interation = zeros(length(intial),3);
for i=1:length(intial)
    x_initial = intial(i);
    res = numerical_sqrt(a,x_initial);
    interation(i,:) = res(2:end);
end
interation(1:5,:)

plot(intial,interation(:,1),'r.') ;  hold on
 plot(intial,interation(:,2),'b.'); hold on
plot(intial,interation(:,3),'g.'); hold on 
xlabel('Initial values x0');
ylabel('n/iteration times');
legend('Newton=x^2-a','Newton=1-a/x^2','Third order');
title('Iteration times with diff initial values sqrt(0.111)');


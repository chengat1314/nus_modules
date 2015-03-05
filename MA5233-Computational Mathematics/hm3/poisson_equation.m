function poisson_equation()
clear,clc
N = 10;
ax = 0;bx = 1;
ay = 0;by = 1;
hx = (bx - ax)/N;
hy = (by - ay)/N;
f_xy = @(x,y)(sin(5*pi*x) * sin(7*pi*y)); %target function 
b = zeros(N^2,1); k = 1;
for xi = ax + hx:hx:bx
    for yi = ay + hx:hy:by
        b(k) = hx*hy*f_xy(xi,yi);
        xx(k) = xi;yy(k) = yi;
        k = k + 1 ;
    end
end 
t = 1; %iteration counter
% steep_method = zeros(N,1);
  gc_method = zeros(N,1);

% for i = 1:N
    A = two_dim_poisson(N);
    % size(A),size(b)
%     sprintf('get the matrix A level at %i agree', i) ;
      [cc,k2] = Conjugate_Gradient(A,b,10000);
%     [ss,k1] = Steepest_decent(A,b,1000000);
%     steep_method(t) = k1;
      gc_method(t) = k2;
    t = t + 1 ;
% end

% plot([10:2:40],steep_method,'r*-');
% hold on
% plot([10:2:40],gc_method,'g^-');

%% plot the result
% zzz = zeros(N,N);
% for i = 1:N
%     for j = 1:N
%         zzz(i,j) = cc((i-1)*N + j) ;
%     end
% end
% figure(1)
% scatter3(xx,yy,cc)
% figure(2)
% 
% [xxx,yyy] = meshgrid(ax + hx:hx:bx);
% contour3(xxx,yyy,zzz,50)
% 

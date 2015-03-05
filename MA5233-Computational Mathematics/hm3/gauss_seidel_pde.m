function gauss_seidel_pde()
%% Programming the Gauss Siedel Method
% Generating the solution matrix
clc,clear
N = 100;
n = N; m= N;
zz=zeros(n+1,m+1);
f=@(x,y)(sin(5*pi*x)*sin(7*pi*y));
tolerance = 10^(-5) ; %tolerance
% Assigning the temperature boundary conditions
h = 1/m;
iteration=0;
count=(m-1)*(n-1);
while(count>0)
count=0;
    for i=2:1:n
        for j=2:1:m
            old=zz(i,j);
            zz(i,j)=( zz(i+1,j)+zz(i-1,j)+zz(i,j+1)+zz(i,j-1) + h^2*f((i-1)*h,(j-1)*h) )/4;
            error(i,j)= abs((zz(i,j)-old)); 
        end 
    end
    if(sum(sum(error))>tolerance)
        count=1;
    end    
iteration=iteration+1; 
end 
iteration
%% Plotting
k = 1;
for i=1:m+1
    for j=1:n+1
        xx(k) = (i-1)*h;
        yy(k) = (j-1)*h;
        zzz(k) = zz(i,j);
        k = k+1;
    end
end 
%% plot colorbar
x1d = 0:h:1;
y1d = 0:h:1;
figure(1)
h_surf=surf(x1d,y1d,zz);view(2)
shading interp
grid on
colorbar   
%% plot 2d point and contour3
figure(2)
scatter3(xx,yy,zzz)

figure(3)
[xxx,yyy] = meshgrid(0:h:1);
contour3(xxx,yyy,zz,50)



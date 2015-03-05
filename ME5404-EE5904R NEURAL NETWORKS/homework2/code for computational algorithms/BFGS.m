x=[-1.2;1.2];
H=diag([2,6]);

tol = 10^(-20);

[y,grad]=rosenbrock(x);
dist=2*tol;

iter=0;
while dist>tol
    
    [value,grad]=rosenbrock(x);
    p=-H*grad;
    
    alpha=linesearch_secant(@rosenbrock,p,x);
    
    iter=iter+1;
    x=x+alpha*p;
      
    s=alpha*p;

    dist=norm(s);
    
    [newvalue,newgrad]=rosenbrock(x);

    y = newgrad-grad;
    
    rho=1/(y'*s);
    
    H=(eye(2)-rho*s*y')*H*(eye(2)-rho*y*s')+rho*s*s';
    
end

iter

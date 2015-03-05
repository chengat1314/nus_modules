function alpha = linesearch_symb(func, d, x)

syms t
f=rosenbrock(x+t*d);
tvalues=solve(diff(f,t));
alpha = real(double(tvalues));
alpha=sort(alpha);
alpha=alpha(find(alpha>0));
alpha=alpha(1);
end
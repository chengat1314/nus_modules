
function fxy = df(x,y) %% 
dfx = 2*(x-1)-400*(y-x^2)*x;
dfy = 200*(y-x^2);
fxy = [dfx,dfy];
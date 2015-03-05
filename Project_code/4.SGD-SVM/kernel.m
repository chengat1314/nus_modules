function k = kernel(x, y);
% function k = kernel(x, y);
%
%	x: (Lx,N) with Lx: number of points; N: dimension
%	y: (Ly,N) with Ly: number of points
%	k: (Lx,Ly)
%
%       KTYPE KSCALE are global and need to be defined
%
%	KTYPE = 1:      linear kernel:      x*y'
%	KTYPE = 2,3,4:  polynomial kernel:  (x*y'*KSCALE+1)^KTYPE
%	KTYPE = 5:      sigmoidal kernel:   tanh(x*y'*KSCALE)
%	KTYPE = 6:	gaussian kernel with variance 1/(2*KSCALE)
%
%       assumes that x and y are in the range [-1:+1]/KSCALE (for KTYPE<6)

global KTYPE
global KSCALE

k = x*y';
if KTYPE == 1				% linear
    % take as is
elseif KTYPE <= 4			% polynomial
    k = (k*KSCALE+1).^KTYPE;
elseif KTYPE == 5			% sigmoidal
    k = tanh(k*KSCALE);
elseif KTYPE == 6			% gaussian
    [Lx,N] = size(x);
    [Ly,N] = size(y);
    k = 2*k;
    k = k-sum(x.^2,2)*ones(1,Ly);
    k = k-ones(Lx,1)*sum(y.^2,2)';
    k = exp(k*KSCALE);
end

function z = cumtrapz_adj(x,y,dim)
%CUMTRAPZ Cumulative trapezoidal numerical integration.
%   Z = CUMTRAPZ(Y) computes an approximation of the cumulative
%   integral of Y via the trapezoidal method (with unit spacing).  To
%   compute the integral for spacing different from one, multiply Z by
%   the spacing increment.
%
%   For vectors, CUMTRAPZ(Y) is a vector containing the cumulative
%   integral of Y. For matrices, CUMTRAPZ(Y) is a matrix the same size as
%   X with the cumulative integral over each column. For N-D arrays,
%   CUMTRAPZ(Y) works along the first non-singleton dimension.
%
%   Z = CUMTRAPZ(X,Y) computes the cumulative integral of Y with respect
%   to X using trapezoidal integration.  X and Y must be vectors of the
%   same length, or X must be a column vector and Y an array whose first
%   non-singleton dimension is length(X).  CUMTRAPZ operates across this
%   dimension.
%
%   Z = CUMTRAPZ(X,Y,DIM) or CUMTRAPZ(Y,DIM) integrates along dimension
%   DIM of Y. The length of X must be the same as size(Y,DIM)).
%
%   Example:
%       Y = [0 1 2; 3 4 5]
%       cumtrapz(Y,1)
%       cumtrapz(Y,2)
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   See also TRAPZ, CUMSUM, INTEGRAL.

%   Copyright 1984-2016 The MathWorks, Inc. 

%   Make sure x and y are column vectors, or y is a matrix.

perm = []; nshifts = 0;
if nargin == 3, % cumtrapz(x,y,dim)
    if ~isscalar(dim) || ~isnumeric(dim)
        error(message('MATLAB:getdimarg:dimensionMustBePositiveInteger'));
    end
    perm = [dim:max(length(size(y)),dim) 1:dim-1];
    y = permute(y,perm);
    [m,n] = size(y);
elseif nargin==2 && isscalar(y) % cumtrapz(y,dim)
    dim = y; y = x;
    if ~isnumeric(dim)
        error(message('MATLAB:getdimarg:dimensionMustBePositiveInteger'));
    end
    perm = [dim:max(length(size(y)),dim) 1:dim-1];
    y = permute(y,perm);
    [m,n] = size(y);
    x = 1:m;
else % cumtrapz(y) and cumtrapz(x, y)
    if nargin < 2
        y = x;
    end
    [y,nshifts] = shiftdim(y);
    [m,n] = size(y);
    dim = nshifts + 1;
    if nargin < 2
        x = 1:m;
    end
end
if ~isvector(x)
    error(message('MATLAB:cumtrapz:xNotVector'));
end
% Make sure we have a column vector.
x = x(:);
if length(x) ~= m
    error(message('MATLAB:cumtrapz:LengthXMismatchY',dim));
end

if isempty(y)
    z = y;
else
    dt = repmat(diff(x,1,1)/2,1,n);
    z = [zeros(1,n,class(y)); cumsum(abs(dt) .* abs((y(1:m-1,:) + y(2:m,:))),1)]; %#ok<ZEROLIKE>
end

siz = size(y);
z = reshape(z,[ones(1,nshifts),siz]);
if ~isempty(perm)
    z = ipermute(z,perm);
end

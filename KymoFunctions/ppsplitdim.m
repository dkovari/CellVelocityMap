function varargout = ppsplitdim(pp)
%break an N-dimensional pp-form spline into N seperate pp splines
% Input:
%   pp: pp-form spline structure
% Output:
%   [pp1,pp2,...] = pp-form splines for each dimension in pp


if nargout > pp.dim
    error('the supplied pp spline only has %d dim, but you specified %d outputs',pp.dim,nargout);
end

DIM = pp.dim;
for d=1:nargout
    p = pp;
    p.coefs = p.coefs(d:DIM:end,:);
    p.dim = 1;
    varargout{d} = p;
end

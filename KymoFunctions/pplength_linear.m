function [cumL,S] = pplength_linear(pp)
%estimate length of a poly spline
%output:
% cumL: cumulative length
% S: points in pp at which cumL is evaluated

npts = 2*numel(pp.breaks);
S = linspace(pp.breaks(1),pp.breaks(end),npts);

X = ppval(pp,S);

dX = diff(X,1,2);

cumL = cumsum(sqrt(sum(dX.^2,1)));
cumL = [0,cumL];
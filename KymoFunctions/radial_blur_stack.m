function res = radial_blur_stack(stack,R,varargin)
% Radially blur an image stack using a specified method/function
% Inputs:
%	stack:	input image stack, single value or RGB both work
%	R:	radius of blur region
% Optional Parameters:
%	'method':
%		'linear' uses a symetric linear kernel k=(R-abs(r))
%		'parabolic' parabolic kernel k=R^2-r^2
%		'hat'	a flat hat function
%		'gaussian' convolve with a gaussian of sigma=R/2
%		@fcn(r): kernel generated using a user defined radial function
%	'shape'
%		'same'	(default) uses zero padding, return same size as I, see conv2()
%		'pad'	uses edge padding, returns same size as I
%		'valid' no padding, returns only central region, see conv2()
%		'full' returns full convolution, see conv2()
% Dependancies:
%   radial_blur() by Daniel T. Kovari
% See radial_blur()for more info.

if ndims(stack)>3
    error('bwstack must be 3d, [H,W,nFrames]')
end

res = zeros(size(stack));

for f = 1:size(stack,3)
    res(:,:,f) = radial_blur(stack(:,:,f),R,varargin{:});
end

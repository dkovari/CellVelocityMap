function data = reinit(data,varargin)%,accuracy, g, ncells)
% reinit: reinitialized data using signedDistanceIterative.
%   data = reinit(initialType, accuracy, grid)
%
% Inputs:
%   data:   Implicit surface function
% Optional Inputs:
%   'accuracy'     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   'grid',g            Grid structure on which data was computed.
%   'NumCells',[d1,d2,...]
%                   number of cell along each dimension to 
%                   reinitialize, defaults to size of grid/data
%==========================================================================
% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
%% Change Log
% Ian Mitchell, 2/14/04
% Dan Kovari 2015-01-31: modified for use with propagateLS...

%% Initialize LS toolbox
% Make sure we can see the kernel m-files.
run('addPathToKernel');

% Across how many grid cells should we reinitialize?
%   (Choose inf to reinitialize to completion).
reinitGridCells = inf;

% What is the convergence criterion?
errorMax = 1e-3;

%% Default parameters.
p = inputParser;
p.CaseSensitive=false;
addParamValue(p,'accuracy','medium',@isaccuracy);
addParamValue(p,'grid',[]);
addParamValue(p,'NumCells',Inf);
parse(p,varargin{:});

accuracy = p.Results.accuracy;
g = p.Results.grid;
if isempty(g)
    %setup default grid
    g.dim = ndims(data);
    g.N = reshape(size(data),[],1);
    g.min = -ones(g.dim,1);
    g.max = ones(g.dim,1);
    g.bdry = @addGhostExtrapolate;
    g = processGrid(g); %process the grid
end

ncells = p.Results.NumCells(1);



%---------------------------------------------------------------------------
% Choose maximum time to reinitialize the right number of grid cells.
%   Of course, there is no point in reinitializing longer than the 
%   diameter of the computational domain.
tMax = min(ncells*max(g.dx),norm(g.max - g.min));

%---------------------------------------------------------------------------
% Reinitialize to signed distance function.
startTime = cputime;
fprintf('reinitalizing: ');
data = signedDistanceIterative(g, data, accuracy, tMax, errorMax);
endTime = cputime;
fprintf('completed in: %g seconds\n', endTime - startTime);


function res = isaccuracy(s)
res=false;
if ischar(s)
    res = any(strcmpi(s,{'low','medium','high','veryHigh'}));
end


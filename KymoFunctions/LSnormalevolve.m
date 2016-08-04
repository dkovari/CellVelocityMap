function [data,T,r] = LSnormalevolve(data0,target,varargin)
% Transform one level-set to another using evolution along the normal to
% the zero level-set.  See chapter 6 of Osher&Fedkiw for mathematical
% derivation.  This code is derived from normalStarDemo(...) in ToolboxLS
% by Ian Mitchell (http://www.cs.ubc.ca/~mitchell/ToolboxLS/index.html).
%
% Required components:
%   ToolboxLS v1.1.1 or higher
%   addPathToKernel.m is properly formatted and in current folder or path
%   reinit.m based on reinitDemo in ToolboxLS
%
% Inputs:
%   data0:  double array specifying data of initial levelset
%   target: double array specifying target levelset
% Parameters:
%   'accuracy': Controllls the order of approximations
%       'low'         Use odeCFL1 and upwindFirstFirst.
%   	'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%   	'high'        Use odeCFL3 and upwindFirstENO3.
%    	'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   'displayType':    Specifying appearance of displayed results
%       'none'      No display (default)
%       See visualizeLevelSet() for more options.
%   'speedFn':  Specify F(x,t) to use in normal velocity Diff.Eq.
%       'difference': Levelset differencing, Machacek et al. EQ 10(default)
%       'asinh':    BiophysJ:2006:90:1439-52, EQ 11
%       'curvature':    BiophysJ:2006:90:1439-52, EQ 12
%       'asinh_curv':   BiophysJ:2006:90:1439-52, hybrid equation (Fig 6E)
%       'tanh': tanh(differece)
%       'tanh08': tanh(0.8*difference)
%       @function:  specify any function of the form
%           out = speedFunction(t,data,schemeData)
%               where schemeData is a structure containing target and grid
%   'residual': minimum residual before stopping calculation, as determined
%               by r = sqrt(sum(sum(D.^2)))/nnz(D); with D=target-data
%               basically a per-pixel error near the levelset
%               default=10^-4
%   'linearize':   linearize gradient data before processing, see O&F
%       'none', 'data', 'target', 'both'
%   'boundary': Boundary method, see ToolboxLS
%   'dX':   grid spacing
%   'minX': grid upper left corner, or dX/2 if not specified
%       NOTE: grid dimensions are set by size(data0);
%   'OutputStack'   if true and ndim < 4, output of intermediate steps will
%                   be saved to data, which has dim=ndim+1, T will be a
%                   vector listing the time points (default: true). Set
%                   this to false if you just want to plot the level set
%                   evolution (also need to pass options for visualizeLevelSet()
%
% Outputs:
%   data:   transformed levelset
%   tNow:   time at exit of function
%
% Copyright Daniel T. Kovari & Wenbin Wei, 2013. All rights reserved.

%-------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('addPathToKernel');

p = inputParser;

addParamValue(p,'accuracy','medium',@isaccuracy);
addParamValue(p,'displayType','none',@isstr);
addParamValue(p,'speedFn','difference',@isspeedFn);
addParamValue(p,'timeStep',0.1,@isnumeric);
addParamValue(p,'t0',0,@isnumeric);
addParamValue(p,'nSteps',20,@isnumeric);
addParamValue(p,'plotSteps',5,@isnumeric);
addParamValue(p,'residual',10^-4,@isnumeric);
addParamValue(p,'linearize','none',@islinearize);
addParamValue(p,'boundary','Extrapolate',@isboundary)
addParamValue(p,'dX',[],@isnumeric);
addParamValue(p,'minX',[],@isnumeric);
addParamValue(p,'maxX',[],@isnumeric);
addParamValue(p,'useSubplots',true,@islogical);
addParamValue(p,'deleteLastPlot',false,@islogical);
addParamValue(p,'OutputStack',true,@islogical);

parse(p,varargin{:});
%set param variables
accuracy = p.Results.accuracy;
displayType = p.Results.displayType;
dt = p.Results.timeStep;
tNow = p.Results.t0;
nSteps = p.Results.nSteps;
plotSteps = p.Results.plotSteps;
residual = p.Results.residual;
T=tNow+dt;

%setup speed function
if isa(p.Results.speedFn,'function_handle')
    speedFn = p.Results.speedFn;
else
    switch(p.Results.speedFn)
        case 'asinh_curv'
            speedFn = @asinh_curv_speedFn;
        case 'difference'
            speedFn = @difference_speedFn;
        case 'asinh'
            speedFn = @asinh_speedFn;
        case'curvature'
            speedFn = @curvature_speedFn;
        case 'tanh'
            speedFn = @tanh_speedFn;
        case 'tanh08'
            speedFn = @tanh08_speedFn;
        case 'tanh06'
            speedFn = @tanh06_speedFn;
        otherwise
            error('unknown speedFn')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup grid
if any(size(data0)~=size(target)) %setup grid dimensions
    error('data0 and target must be same dimension');
end
g.dim = ndims(data0);
g.N = reshape(size(data0),[],1);
if ~isempty(p.Results.dX)   %setup grid resolution
    if (ndims(p.Results.dX)>1&&ndims(p.Results.dX)~=ndims(data0))
        error('dX must be either scalar or same ndims as data0');
    end
    g.dx = reshape(p.Results.dX,[],1);
end

if ~isempty(p.Results.minX) %setup grid min/max values
    if ndims(p.Results.minX)>1&&ndims(p.Results.minX)~=ndims(data0)
        error('minX must be scalar or same ndims as data0');
    end
    g.min = reshape(p.Results.minX,[],1);
    if ~isempty(p.Results.dX)
        g.max = g.min+(g.N-1).*g.dx;
    else
        if isempty(p.Results.maxX)
            error('if minX is specified maxX or dX must be specified');
        else
            g.max = reshape(p.Results.maxX,[],1);
        end
    end
else
    if ~isempty(p.Results.dX)
        g.min = ones(g.dim,1)*g.dx/2;
        g.max = g.min+(g.N-1).*g.dx;
    else
        g.min = -ones(g.dim,1);
        g.max = ones(g.dim,1);
    end
end

switch(p.Results.boundary)  %setup boundary conditions
    case 'Extrapolate'
        g.bdry = @addGhostExtrapolate;
    case 'Periodic'
        g.bdry = @addGhostPeriodic;
    otherwise
        error('boundary type not understood');
end

g = processGrid(g); %process the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------
% Setup Stack


if(p.Results.OutputStack)
    OutputStack = true;
    switch(g.dim)
        case 1
            stack = nan(numel(data0),nSteps);
            T=nan(nSteps,1);
        case 2
           stack = nan([size(data0),nSteps]);
           T=nan(nSteps,1);
        case 3
            stack = nan([size(data0),nSteps]);
            T=nan(nSteps,1);
        otherwise
            OutputStack = false;
    end
else
    OutputStack = false;
end
%------------------------------------------------------------------------


%reinitialize data
switch(p.Results.linearize)
    case 'target'
        data=data0;
        target=reinit(target,'accuracy','medium','grid',g);
    case 'data'
        data=reinit(data0,'accuracy','medium','grid',g);
    case 'both'
        target=reinit(target,'accuracy','medium','grid',g);
        data=reinit(data0,'accuracy','medium','grid',g);
    case 'none'
        data = data0;
    otherwise
        error('reinit not understood')
end
%--------------------------------------------------------------------------
% Set up motion in the normal direction (derivative choice is set below).
% Choose approximations at appropriate level of accuracy.
switch(accuracy)
    case 'low'
        schemeData.derivFunc = @upwindFirstFirst;
        integratorFunc = @odeCFL1;
    case 'medium'
        schemeData.derivFunc = @upwindFirstENO2;
        integratorFunc = @odeCFL2;
    case 'high'
        schemeData.derivFunc = @upwindFirstENO3;
        integratorFunc = @odeCFL3;
    case 'veryHigh'
        schemeData.derivFunc = @upwindFirstWENO5;
        integratorFunc = @odeCFL3;
    otherwise
        error('Unknown accuracy level %s', accuracy);
end

schemeFunc = @termNormal;
schemeData.grid = g;
schemeData.speed = speedFn;
schemeData.target=target;
schemeData.terminal=5;
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');
if(nSteps<2)
	integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%--------------------------------------------------------------------------
% Initialize Display
if ~strcmpi('none',displayType)
    level=0;
    figure();
    h_data0 = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    if(g.dim > 1)
        axis(g.axis);
        axis equal;
    end
    title('Initial Data');
    
    figure();
    h_target = visualizeLevelSet(g, target, displayType, level, [ 't = ' num2str(tNow) ]);
    if(g.dim > 1)
        axis(g.axis);
        axis equal;
    end
    title('Target Data');
    
    f_data=figure();
    if p.Results.useSubplots
        rows = ceil(sqrt(plotSteps));
        cols = ceil(plotSteps / rows);
        plotNum = 1;
        subplot(rows, cols, plotNum);
    end
    h_data = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    if(g.dim > 1)
        axis(g.axis);
        axis equal;
    end
    hold on;
end
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Loop for nSteps or until residual is reached
%setup residual
%only care about region near zero levelset
%assume levelset is approx a signed dist with slope 1
%preserve about +/- 5 pixels on each side
rw= 5*sqrt(reshape(g.dx,1,[])*reshape(g.dx,[],1));
r=Inf;
rTar = target;
rTar(rTar>rw)=rw;
rTar(rTar<-rw)=-rw;
% B = rTar>-rw&rTar<rw;
% figure(99);
% imagesc(B);
nPx = nnz(rTar>-rw&rTar<rw);
for istep = 1:nSteps
    %fprintf('step: %d\n',istep);
    %calculate residual
    if r<residual
        break;
    end
    % Take a timestep.
    [t,y] = feval(integratorFunc, schemeFunc,...
        [tNow,tNow+dt], data(:),...
        integratorOptions, schemeData);
    data=reshape(y,g.shape);    % Get back the correctly shaped data array
    tNow = t(end);  %set current computation time
    
    if OutputStack
        T(istep)=tNow;
        switch(g.dim)
            case 1
                stack(:,istep)=data(:);
            case 2
                stack(:,:,istep)=data;
            case 3
                stack(:,:,:,istep)=data;
        end
    end
    
    %plot if necessary
    if ~strcmpi('none',displayType)&& mod(istep,nSteps/plotSteps)==0
        figure(f_data); %get figure;
        [ figure_az, figure_el ] = view; %store previous view

        if(p.Results.deleteLastPlot)  % Delete last visualization if necessary.
            delete(h_data);
        end

        if(p.Results.useSubplots) % Move to next subplot if necessary.
            plotNum = plotNum + 1;
            subplot(rows, cols, plotNum);
        end
        % Create new visualization.
        h_data = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
        % Restore view.
        view(figure_az, figure_el);
    end
    
    %calculate residual
    rDat =data;
    rDat(data>rw)=rw;
    rDat(data<-rw)=-rw;
    D = rDat-rTar;
%     numel(D)
%     nnz(D)
%     nPx
%     figure(100);clf;
%     imagesc(D.^2);
%     sqrt(sum(sum(D.^2)))
    %r = sqrt(trace(D'*D));
    r = sqrt(sum(sum(D.^2)))/nPx;
    fprintf('residual: %f\n',r);
end
if ~strcmpi('none',displayType)
    figure();
    visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    if(g.dim > 1)
        axis(g.axis);
        axis equal;
    end
    title('Final Data');
end

if OutputStack
    T(istep:end)=[];
    switch(g.dim)
        case 1
            stack(:,istep:end)=[];
        case 2
            stack(:,:,istep:end)=[];
        case 3
            stack(:,:,:,istep:end)=[];
    end
    data=stack;
else
    T=tNow;
end


function out = asinh_curv_speedFn(~,data,schemeData)
compFact = 1;
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
[curvature,~] = curvatureSecond(schemeData.grid,data);
out=((data-target)>=0).*asinh(data-target)+((data-target)<0).*(data-target).*(compFact+curvature);

function out = asinh_speedFn(~,data,schemeData)
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
out=asinh(data-target);

function out = difference_speedFn(~,data,schemeData)
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
out=data-target;

function out = tanh_speedFn(~,data,schemeData)
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
out=tanh(data-target);

function out = tanh08_speedFn(~,data,schemeData)
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
out=tanh(0.8*(data-target));

function out = tanh06_speedFn(~,data,schemeData)
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
out=tanh(0.6*(data-target));

function out = curvature_speedFn(~,data,schemeData)
compFact = 1/3;
checkStructureFields(schemeData, 'target','grid');
target=schemeData.target;
[curvature,~] = curvatureSecond(schemeData.grid,data);
out=(data-target).*(compFact+curvature);

function res = isspeedFn(s)
res = false;
if isa(s,'function_handle')
    res = true;
    return
end
if ischar(s)
    res = any(strcmpi(s,{'asinh_curv','difference','asinh','curvature','tanh','tanh08','tanh06'}));
end

function res = isaccuracy(s)
res=false;
if ischar(s)
    res = any(strcmpi(s,{'low','medium','high','veryHigh'}));
end

function res = islinearize(s)
res=false;
if ischar(s)
    res = any(strcmpi(s,{'none','both','data','target'}));
end

function res = isboundary(s)
res=false;
if ischar(s)
    res = any(strcmpi(s,{'Extrapolate','Periodic'}));
end
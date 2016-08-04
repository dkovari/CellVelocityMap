function [pp,range,C] = implicit2explicit(data,varargin)
%% Function Description
% Converts the implicit function specified by the datapoints in data to an
% explicit parametric function specified by the zero crossing of that
% function.
%
% Requirements:
%   contourcs() by Takeshi Ikuma
%   http://www.mathworks.com/matlabcentral/fileexchange/28447-contourcs-to-obtain-contourc-output-as-a-struct-array
%
% Input
%   data: A matrix [H,W] specifying datapoints to use for the surface
% Parameters
%   'ForcePeriodic' Logical specifying if parametric curve should be
%                   periodic over the range of its input. (default true)
%   'LevelSet'  value specifying which levelset contour to use (def. 0)
%   'dX'    x spacing of gridpoints (default 1)
%   'dY'    y spacing of gridpoints (default 1)
%   'Xmin'  initial x value of grid (default 1*dX)
%   'Ymin'  initial y value of gris (default 1*dX)
%   'Xmax'  end value of x grid (default = size(data,2) or set by dX&Xmin)
%   'Ymax'  end value of y grid (default = size(data,1) or set by dY&Ymin)
%   'Xpoints'	vector specifying x locations
%   'Ypoints'	vector specifying y locations
%   'Xgrid'     matrix specifying xgrid
%               Note: contourc only accepts a vector specifying the grid
%               dims so only the first vector Xgrid(1,:) is used
%   'Ygrid'     matrix specifying ygrid
%               Note: contourc only accepts a vector specifying the grid
%               dims so only the first vector Ygrid(:,1) is used
%   'LargestOnly', t/f  Return only the spline for the longest contour
%   'SortLength', 'ascend', 'descend', 'none'
%                   Sort the splines by length, default: none (unsorted)
%
% Outputs
%   pp:     PP-form of the spline-fit explicit parametric function
%   range:  [min(pp.breaks),max(pp.breaks)] specifying range or function
%   C: array of points on the level set created by contourcs()
% Note:
%   If multiple contours are found, pp is a column vector and range is a
%   matrix with row cooresponding to rows in pp.  They are aranged in
%
% Copyright Daniel T. Kovari & Wenbin Wei, 2013. All rights reserved.

%% Parse Input options
% 
p=inputParser;
p.KeepUnmatched=true;
p.CaseSensitive=false;
addParamValue(p,'Xgrid',[],@isnumeric);
addParamValue(p,'Ygrid',[],@isnumeric);
addParamValue(p,'dX',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'dY',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'Xmin',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'Xmax',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'Ymin',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'Ymax',NaN,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'Xpoints',[],@isnumeric);
addParamValue(p,'Ypoints',[],@isnumeric);
addParamValue(p,'ForcePeriodic',true,@islogical);
addParamValue(p,'LevelSet',0,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'LargestOnly',false,@islogical);
addParamValue(p,'SortLength','none',@(x) ischar(x)&&(strcmpi(x,'ascend')||strcmpi(x,'descend')));
parse(p,varargin{:});

if ~ismatrix(data)
    error('only works with 2D data');
end
x=1:size(data,2);
y=1:size(data,1);
dX=1;
dY=1;
if ~isnan(p.Results.dX)
    dX=p.Results.dX;
    x=x*dX;
end
if ~isnan(p.Results.dY)
    dY=p.Results.dY;
    y=y*dY;
end
if ~isnan(p.Results.Xmin)
    x=x-x(1)+p.Results.Xmin;
end
if ~isnan(p.Results.Ymin)
    y=y-y(1)+p.Results.Ymin;
end
if ~isnan(p.Results.Xmax)
    if ~isnan(p.Results.Xmin)
        x=linspace(x(1),p.Results.Xmax,size(data,2));
    else
        Xmax = p.Results.Xmax;
        x=linspace(Xmax-sixe(data,2)*dX,Xmax,size(data,2));
    end
end
if ~isnan(p.Results.Ymax)
    if ~isnan(p.Results.Ymin)
        y=linspace(y(1),p.Results.Ymax,size(data,1));
    else
        Ymax = p.Results.Ymax;
        y=linspace(Ymax-sixe(data,1)*dY,Ymax,size(data,1));
    end
end
if ~isempty(p.Results.Xpoints)
    if numel(p.Results.Xpoints)~=size(data,2)
        error('numel(Xpoints) must be same as width of data');
    end
    x=p.Results.Xpoints;
end
if ~isempty(p.Results.Ypoints)
    if numel(p.Results.Ypoints)~=size(data,1)
        error('numel(Ypoints) must be same as height of data');
    end
    y=p.Results.Ypoints;
end
%[XX,YY]=meshgrid(x,y);
if ~isempty(p.Results.Xgrid)
    if ~all(size(p.Results.Xgrid)==size(data))
        error('dim of Xgrid must match data');
    end
    XX=p.Results.Xgrid;
    x=XX(1,:);
end
if ~isempty(p.Results.Ygrid)
    if ~all(size(p.Results.Ygrid)==size(data))
        error('dim of Ygrid must match data');
    end
    YY=p.Results.Ygrid;
    y=YY(:,1);
end
vv = ones(1,2)*p.Results.LevelSet;
ForcePeriodic = p.Results.ForcePeriodic;
LargestOnly = p.Results.LargestOnly;
SortLength = p.Results.SortLength;
%End Input Parsing
%% Process Data
% DESCRIPTIVE TEXT
C=contourcs(x,y,data,vv);   %get contours from image, output in sturt form
%note: C needs to be a column vector, in current version of contourcs it is
if(numel(C)<1)
    error('No contour found')
end

if LargestOnly
    L = [C(:).Length];
    [~,idx] = max(L);
    C=C(idx);
else
    if ~strcmpi(SortLength,'none')
        L = [C(:).Length];
        [~,idx] = sort(L,SortLength);
        C = C(idx);
    end
end
if ForcePeriodic
    C = arrayfun(@MakePeriodic,C);
end
[pp,range]=arrayfun(@splinefitc,C,'UniformOutput',false);
pp = cell2mat(pp);
range = cell2mat(range);

function C = MakePeriodic(C)
if C.X(end)~=C.X(1)||C.Y(end)~=C.Y(1)
    C.Length = C.Length+1;
    C.X(end+1)=C.X(1);
    C.Y(end+1)=C.Y(1);
end
    
function [pp,range] = splinefitc(C)
pp = cscvn([reshape(C.X,1,[]);reshape(C.Y,1,[])]);
range=[min(pp.breaks),max(pp.breaks)];

        

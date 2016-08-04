function [S,Xi,Yi] = intersectSpline2d(pp,slim,X,Y,varargin)
% Finds the exact intersection of a line series <X,Y> with a spline
% function.
% 
% Requirements:
%   intersections() by Douglas Schwarz 
%   http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
%
% Inputs
%   pp: PP-form spline function
%   slim: limits overwhich to search Spline
%   X: vector of X coordinates for second line
%   Y: vector of Y coordinates for seconf line
% Parameters
%   'ppResolution'  default 1000, initial resolution of spline
%   'ForcePeriodic' default true, forces spline to be interpreted as
%                   periodic over its range min,max(pp.breaks)
%
% Copyright Daniel T. Kovari & Wenbin Wei, 2013. All rights reserved.


p = inputParser;
p.KeepUnmatched=true;
p.CaseSensitive=false;

addParamValue(p,'ppResolution',1000,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'ForcePeriodic',true,@islogical);


parse(p,varargin{:});
ppRes = p.Results.ppResolution;
periodic = p.Results.ForcePeriodic;

%evaluate pp over ppRes points between slim
SS = linspace(slim(1),slim(2),ppRes);

ppXY = fnval(pp,SS);

%find approximate
[Xi,Yi]=intersections(X,Y,ppXY(1,:),ppXY(2,:));

if isempty(Xi)
    S = [];
    return
end
S=nan(size(Xi));
for f = 1:numel(Xi)
    %find points in ppXY closest to <Xi,Yi>
    [~,idx]=min(sqrt( (ppXY(1,:)-Xi(f)).^2 + (ppXY(2,:)-Yi(f)).^2));
    S(f)=SS(idx);
    
    %find exact intersection
    if numel(X)>2 %find points on line closest to approx intersection
        [~,Lidx1] = min( sqrt( (Xi(f)-X(:)).^2 + (Yi(f)-Y(:)).^2 ) );
        LX1 = X(Lidx1);
        LY1 = Y(Lidx1);
        tmpX=X;
        tmpY=Y;
        tmpX(Lidx1)=[];
        tmpY(Lidx1)=[];
        [~,Lidx2] = min( sqrt( (Xi(f)-tmpX(:)).^2 + (Yi(f)-tmpY(:)).^2 ) );
        LX2 = tmpX(Lidx2);
        LY2 = tmpY(Lidx2);
    else
        LX1=X(1);
        LY1=Y(1);
        LX2=X(2);
        LY2=Y(2);
    end
    
    fn = @(t) distfun(t,pp,[LX1;LX2],[LY1;LY2],periodic);
    
    if periodic
        thisVal = fn(SS(idx));
        if thisVal~=0
            switch idx
                case 1
                    ds = SS(idx+1)-SS(idx);
                    ns = SS(idx+1);
                    ps = SS(idx)-ds;
                case length(SS)
                    ds = SS(idx)-SS(idx-1);
                    ns = SS(idx)+ds;
                    ps = SS(idx-1);
                otherwise
                    ns = SS(idx+1);
                    ps = SS(idx-1);
            end
            if sign(thisVal)~=sign(fn(ns)) && sign(thisVal)~=sign(fn(ps))
                %sign change in both directions therefore guess the closer
                %point
                pxy = fnval(pp,ps);
                nxy = fnval(pp,ns);
                pdis = sqrt( (Xi(f)-pxy(1)).^2 + (Yi(f)-pxy(2)).^2 );
                ndis = sqrt( (Xi(f)-nxy(1)).^2 + (Yi(f)-nxy(2)).^2 );
                if pdis<ndis
                    Srange = [SS(idx), ps];
                else
                    Srange = [SS(idx),ns];
                end
            elseif sign(thisVal)~=sign(fn(ns))
                Srange = [SS(idx),ns];
            elseif sign(thisVal)~=sign(fn(ps))
                Srange = [SS(idx),ps];
            else
                Srange = SS(idx);
            end
                    
        else
            Srange = SS(idx);
        end
    else
        Srange = SS(idx);
    end
    S(f) = fzero(fn,Srange);
    
    xy=fnval(pp,S(f));
    Xi(f)=xy(1);
    Yi(f)=xy(2);
    
end

function d = distfun(t,pp,X,Y,periodic)
if periodic
    t=mod(t,max(pp.breaks));
end
rhat = [X(2)-X(1);Y(2)-Y(1)]/sqrt((X(2)-X(1))^2+(Y(2)-Y(1))^2);
vhat = [-rhat(2);rhat(1)];
XY1 = [X(1);Y(1)];
d = vhat'*(fnval(pp,t(:))-XY1);

        
        


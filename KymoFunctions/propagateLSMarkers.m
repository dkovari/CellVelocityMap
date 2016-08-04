function [markers, PP, SS] =  propagateLSMarkers(stack,varargin)
%%Propegate markers along the outline of BW regions that evolve with time
% Requirements:
%   intersectSpline2d() by Daniel T. Kovari
%   implicit2explicit() by Daniel T. Kovari
%   intersections() by Douglas Schwarz
% Inputs
%   stack:    logical(Y,X,t) 3D array specifying evolving implicit function
% Parameters
%   'nMarkers'  number of markers to propegate, equally spaced along the
%               levelset of stack(:,:,1)
%   'initMarkers' structure specifying initial markers
%               initMarkers.PP = ppform spline specifying shape
%               initMarkers.SS = vector of value overwhich to evaluate PP
%   'DisplayDebug'  Generates animation illustrating process, for debug.
%   See implicit2explicit() for more parameters.
%
% Output
%   markers:    [2,nMarkers,t] array specifying XY by nMarkers points found
%               along the edges markers(1,:,t)->X vals, markers(2,:,t)->Y
%
%   PP:         struct array specifying ppform splines of the levelsets
%               each element in PP is the spline corresponding to the
%               calculated 
%
%   SS:         [nMarkers,t] array specifying marker points along PP{t}
%               NOTE: markers(:,:,t)=fnval(PP(t),SS(:,t));
%
% Copyright Daniel T. Kovari & Wenbin Wei, 2013. All rights reserved.

%% Parse Input options
p=inputParser;
p.KeepUnmatched=true;
p.CaseSensitive=false;
%Parameter Definitions
addParamValue(p,'nMarkers',50,@(x) isnumeric(x)&&isscalar(x));
addParamValue(p,'initMarkers',[],@(x) isstruct(x)||isempty(x));
addParamValue(p,'DisplayDebug',false,@islogical);
addParamValue(p,'LevelSet',0,@(x) isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
LS = p.Results.LevelSet;
DispDebug = p.Results.DisplayDebug;

nF = size(stack,3);
%PP=cell(size(stack,3),1);
PP(nF) = struct('form','pp','breaks',[],'coefs',[],'pieces',0,'order',0,'dim',0);

initMarkers=p.Results.initMarkers;
if isstruct(initMarkers) %use specified markers
    if ~isfield(initMarkers,'SS')||~isfield(initMarkers,'PP')
        error('need field SS & PP in initMarkers struct');
    end
    nMarkers = numel(initMarkers.SS);
    markers = nan(2,nMarkers,size(stack,3));
    SS = nan(nMarkers,size(stack,3));
    prevPP = initMarkers.PP;
    prevSS = reshape(initMarkers.SS,1,[]);
    startF = 1;
    nextRange = [min(prevPP.breaks),max(prevPP.breaks)];
else %determine initial markers
    nMarkers = p.Results.nMarkers;
    markers = nan(2,nMarkers,size(stack,3));
    SS = nan(nMarkers,size(stack,3));
    %find first set of markers
    [prevPP,nextRange]=implicit2explicit(stack(:,:,1),varargin{:});
    prevSS = linspace(nextRange(1),nextRange(2),nMarkers+1);
    prevSS(end)=[];
    %save to first frame
    markers(:,:,1) = fnval(prevPP,prevSS);
    SS(:,1)=prevSS;
    PP(1)=prevPP;
    startF = 2;
end
%End Input Parsing

%% Start Code Here
if DispDebug
    %plot starting LS in green
    dbfig = figure(1); clf; hold on; axis([0,size(stack,2)+1,0,size(stack,1)+1]);
    tmpXY = fnval(prevPP,prevPP.breaks);
    plot(tmpXY(1,:),tmpXY(2,:),'-c');
    %initialize graphics handles, (delete will complain if we dont)
    hsel = [];
    hpre = [];
    hNorm = [];
    hcros = [];
    hcros2 = [];
end
for f = startF:size(stack,3)
    fprintf('propagating to %d\n',f);
    
    %calc normal line
    xy0 = fnval(prevPP,prevSS)'; %x-y position of initial marker
    dxy = fnval(fnder(prevPP),prevSS)';
    dxy_perp = [-dxy(:,2),dxy(:,1)];
    %get intersection with next levelset
    [Xi,Yi,~,PP(f),Si] = intersectLineLevelSet(stack(:,:,f),xy0,dxy_perp,'LevelSet',LS);
    %find closest intersection to previous marker
    for n=1:size(Si,1)
        if isempty(Si{n})
            %we didn't find an intersection with the normal line, just
            %choose closest break point on next curve...this could get
            %buggy
            xy = ppval(PP(f),PP(f).breaks);
            dist = sqrt( (xy(1,:)-xy0(n,1)).^2 + (xy(2,:)-xy0(n,2)).^2 );
            [~,idx] = min(dist);
            SS(n,f) = PP(f).breaks(idx);
            markers(1,n,f) = xy(1,idx);
            markers(2,n,f) = xy(2,idx);
        else
            dist = sqrt( (Xi{n}-xy0(n,1)).^2 + (Yi{n}-xy0(n,2)).^2 );
            [~,idx] = min(dist);
            SS(n,f) = Si{n}(idx);
            markers(1,n,f) = Xi{n}(idx);
            markers(2,n,f) = Yi{n}(idx);
        end
    end
    
    
%     dist = sqrt( (Xi-xy0(:,1)*ones(1,size(Xi,2))).^2 + (Yi-xy0(:,2)*ones(1,size(Yi,2))).^2 );
%     [~,I]=min( dist,[],2 );
%     
%     for n=1:size(Si,1)
%         SS(n,f) = Si(n,I(n));
%         markers(1,n,f) = Xi(n,I(n));
%         markers(2,n,f) = Yi(n,I(n));
%     end
    
    if DispDebug %Debug Display
        figure(dbfig);
        clf;
        imagesc(stack(:,:,f));colormap gray;
        axis image;
        axis xy;
        hold on;
        tmpXY = fnval(PP{f},PP{f}.breaks);
        plot(tmpXY(1,:),tmpXY(2,:),'-b');
        
        hsel = [];
        hpre = [];
        hNorm = [];
        hcros = [];
        hcros2 = [];
        
        for n=1:size(Si,1)
            figure(dbfig);
            delete(hpre);
            delete(hNorm);
            delete(hcros2);
            delete(hcros);
            hpre = plot(xy0(n,1),xy0(n,2),'.g');
            if dxy(n,2)==0    %vertical line spanning LS image
                Xn=[1;1]*xy0(n,1);
                Yn=[0;size(stack,1)+1];
            else
                m = -dxy(n,1)/dxy(n,2);
                Xn = [0;size(stack,2)+1];
                Yn = m*(Xn-xy0(n,1))+xy0(n,2);
                %set line limits to boundary of the LS image
                [Xn,Yn]=intersections(Xn,Yn,[0,size(stack,2)+1,size(stack,2)+1,0,0],[0,0,size(stack,1)+1,size(stack,1)+1,0]);
            end

            hNorm = plot(Xn,Yn,'-r');

            tmpXY = fnval(PP{f},Si(n,:));
            hcros = plot(tmpXY(1,:),tmpXY(2,:),'xr');
            hcros2 = plot(Xi(n,:),Yi(n,:),'om');
            delete(hsel);
            hsel = plot(markers(1,n,f),markers(2,n,f),'sg'); drawnow;
            %pause;
        end
    end
    
    %set "prev" vars for next loop
    prevSS = SS(:,f);
    prevPP = PP(f);
    
end
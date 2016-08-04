function KymoMontage(CellData)
%add sub-functions to path
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
addpath(genpath(fullfile(pathstr,'KymoFunctions')));

TIMES = [2.5, 3.5, 4.5, 5.5]*60;

persistent last_dir;
if nargin<1
%% Load Data
    %select file
    [File,Dir] = uigetfile(fullfile(last_dir,'*.mat'),'Select cell data file');
    if File==0
        return
    end
    if ~isempty(Dir)
        last_dir = Dir;
    end

    CellData = load(fullfile(Dir,File));

elseif ischar(CellData)
    if exist(CellData,'file');
        CellData = load(CellData);
    else
        error('File: %s does not exist',CellData);
    end
elseif isstruct(CellData)
    if nargin<2
        File = '';
        Dir = '';
    else
        [Dir,File,ext] = fileparts(FilePath);
        File = [File,ext];
    end
elseif ~isstruct(CellData)
    error('wrong argument type');
end

if ~any(isfield(CellData,{...
                            'kymo_Lpts',...
                            'kymo_vel',...
                            'kymo_markersF',...
                            'kymo_data',...
                            'kymo_ppF',...
                            'kymo_ssF'}))
    error('Structure doesnt contain kymo variables');
end

%% Get needed data
%Lpts = CellData.kymo_Lpts;
vel = CellData.kymo_vel;
%Time = CellData.Time;
%Area = CellData.Area;
%PxScale = CellData.PxScale;
%L = CellData.kymo_L;
data = CellData.kymo_data;
origstack = CellData.origstack;
ppF = CellData.kymo_ppF;
ssF = CellData.kymo_ssF;
markersF = CellData.kymo_markersF;

if ~isfield(CellData,'PxScale')
    CellData.PxScale = 0.157825;
    CellData.PxScaleUnit = 'µm';
end

if ~isfield(CellData,'PxScaleUnit')
    CellData.PxScaleUnit = 'µm';
end
PX_SCALE = CellData.PxScale;

%% Time Stamp and Scale Bar
%Scale bar size in px
SB_LENGTH = 20; %µm
SB_WIDTH = 4; %figure points
PX_LENGTH = SB_LENGTH/PX_SCALE;
SB_X = [10,10 + PX_LENGTH];
SB_Y = [10,10];

%location of timestamp text
TIMESTAMP_FONT_SIZE = 14;
ylocF = 0.95; %top of text (frac of axis)
xlocF = 0.01; %left of text (frac of axis)

%% Load Colormaps
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
load(fullfile(pathstr,'KymoFunctions','SymetricVelcolormap.mat'),'mycmap')
load(fullfile(pathstr,'KymoFunctions','kymocmap.mat'),'kymocmap');
%cmap = mycmap;


%% Image Montage with arrows
VSCALE = 10/(0.15); %px/(µm/s)
clim=stackclim(origstack,'average');
cmap = gray(255);
H = size(origstack,1);
W = size(origstack,2);

for f=1:numel(TIMES)
    [~,FRAMES(f)] = min(abs(CellData.Time-TIMES(f)));
end

NF = numel(FRAMES);

nPC = ceil(sqrt(NF));
nPR = ceil(NF/nPC);

%% setup fig and axes
hfig = figure('Position',[0,0,nPC*W+20*(nPC-1),nPR*H+20*(nPC-1)]);
set(hfig,'units','inches');
fig_sz = get(hfig,'position');
set(hfig,'papersize',fig_sz(3:4));

%% Colorbar
%Setup color lims for velocity values
climVEL = [-0.15,0.15]; %color limits µm/sec velLim ;

%colorbar ticks
nCTicks = 7;
CBticks = linspace(climVEL(1) ,climVEL(2) ,nCTicks);




for j=1:NF
    f = FRAMES(j);
    hAx = subplot(nPR,nPC,j);
    imbase = ind2rgb( gray2ind( mat2gray(CellData.origstack(:,:,f),clim),size(cmap,1)),cmap);
    image('Parent',hAx,'CData',imbase);
    axis(hAx,'image','off','tight');
    
    
   
    hold on;
    %plot scale bar
    plot(SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);
    
    %% Time stamp
    YLIM = get(hAx,'ylim');
    XLIM = get(hAx,'xlim');

    switch get(hAx,'xdir')
        case 'normal'
            xloc = XLIM(1)+xlocF*(XLIM(2)-XLIM(1));
        case 'reverse'
            xloc = XLIM(2)-xlocF*(XLIM(2)-XLIM(1));
    end
    switch get(hAx,'ydir')
        case 'normal'
            yloc = YLIM(1)+ylocF*(YLIM(2)-YLIM(1));
        case 'reverse'
            yloc = yLIM(2)-ylocF*(YLIM(2)-YLIM(1));
    end

    str = sprintf('Time: %03.01f min',CellData.Time(f)/60);
    text(xloc,yloc,str,'Color','w','FontSize',TIMESTAMP_FONT_SIZE);
    
    set(hAx,'CLim',climVEL);
    %set(hax,'CLimMode','manual');
    %set colormap
    colormap(hAx,mycmap);

    %show colorbar
    hcb = colorbar(hAx,'location','Eastoutside');
    set(hcb,'Ticks',CBticks);

    %set colorbar font
    set(hcb,'FontWeight','bold');
    %colorbar title
    hcbtitle = get(hcb,'ylabel');
    set(hcbtitle,'string','Normal Velocity [µm/s]');
    set(hcbtitle,'FontSize',12);
    
    %% Plot KymoData
    if ~isempty(vel{f})
        %plot perimeter line
        s=linspace(ppF(f).breaks(1),ppF(f).breaks(end),500);
        xy = ppval(ppF(f),s);

        plot(hAx,xy(1,:),xy(2,:),'-','Color','w','LineWidth',1.5);

        %Use gradient to determine inward normal direction
        [Gx,Gy] = gradient(data(:,:,f));
        Vxg = interp2(Gx,markersF{f}(1,1)',markersF{f}(2,1)');
        Vyg = interp2(Gy,markersF{f}(1,1)',markersF{f}(2,1)');


        %calp perp line
        dxy = fnval(fnder(ppF(f),1),ssF{f})';
        Vx = dxy(:,2);
        Vy = -dxy(:,1);

        %do a dot-product to see if vector is outward or inward
        d = dot([Vxg,Vyg],[Vx(1),Vy(1)]);
        if d>0
            Vx = -Vx;
            Vy = -Vy;
        end

        %normalize vectors, make outward normal (V=-V)
        Vx = Vx./sqrt(Vx.^2+Vy.^2);
        Vy = Vy./sqrt(Vx.^2+Vy.^2);

        %scale the velocity
        Vx = Vx.*vel{f}*VSCALE ;
        Vy = Vy.*vel{f}*VSCALE ;

        %plot the arrows in color cooresponding to clim and current colormap
        colors = num2climcolor(mycmap,climVEL,vel{f});
        for k=1:numel(vel{f})
            quiver(hAx,markersF{f}(1,k)',markersF{f}(2,k)',Vx(k),Vy(k),0,'Color',colors(k,:));
        end
    end

end


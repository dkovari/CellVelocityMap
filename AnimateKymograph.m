function AnimateKymograph(CellData,FilePath)

%add sub-functions to path
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
addpath(genpath(fullfile(pathstr,'KymoFunctions')));

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

%% Plot Time Stamp
SHOW_TIMESTAMP = true;
%Time Stamp Placement
TIME_FONT_SIZE = 28;
yloc = 0.99; %top of text (frac of axis)
xloc = 0.01; %left of text (frac of axis)

%% location and size of scalebar
%Scale bar size in px
SB_FONT_SIZE = 28;
SB_LENGTH = 20; %µm
SB_WIDTH = 6; %figure points
PX_LENGTH = SB_LENGTH/(PX_SCALE);

RELATIVE_SB  = true; %place scalebar relative to axes corner
SB_POS = [0.05,0.05]; %position of scalebar [x,y]
SB_X = [10,10 + PX_LENGTH]; %non-relative position of SB
SB_Y = [30,30];%non-relative position of SB

%% Load Colormaps
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
load(fullfile(pathstr,'KymoFunctions','SymetricVelcolormap.mat'),'mycmap')
load(fullfile(pathstr,'KymoFunctions','kymocmap.mat'),'kymocmap');
%cmap = mycmap;

%% Animate with arrows
VSCALE = 10/(0.15); %px/(µm/s)
clim=stackclim(origstack,'average');
cmap = gray(255);
H = size(origstack,1);
W = size(origstack,2);
nF = size(vel,1);

%setup fig and axes
hfig = figure('Position',[0,0,2.5*W+20,2.5*H+20]);%,'Renderer','zbuffer');
hax = axes('Parent',hfig);

%setup image
imbase = ind2rgb( gray2ind( mat2gray(origstack(:,:,1),clim),size(cmap,1)),cmap);
him = image('Parent',hax,'CData',imbase);
set(him,'handlevisibility','callback');
axis(hax,'xy','image','off');


%Setup color lims for velocity values
climVEL = [-0.15,0.15]; %color limits µm/sec velLim ;
set(hax,'CLim',climVEL);
%set(hax,'CLimMode','manual');
%set colormap
colormap(hax,mycmap);

%show colorbar
hcb = colorbar(hax,'location','East');
%colorbar ticks
nCTicks = 7;
CBticks = linspace(climVEL(1) ,climVEL(2) ,nCTicks);
set(hcb,'Ticks',CBticks);
%set(hcb,'YLimMode','manual');
%set(hcb,'YLim',[1,size(mycmap,1)]);
%set(hcb,'YTickMode','manual');
%set(hcb,'YTick',linspace(1,size(mycmap,1),nCTicks));
%set(hcb,'YTickLabel',num2str(CBticks','%0.2f'));

%set colorbar font
set(hcb,'FontWeight','bold');
set(hcb,'YColor','w');
set(hcb,'XColor','w');
%colorbar title
hcbtitle = get(hcb,'ylabel');
set(hcbtitle,'string','Normal Velocity [µm/s]');
set(hcbtitle,'FontSize',12);
%set(hcbtitle,'FontWeight','Bold');

%set(hax,'NextPlot','replacechildren');

if exist('Anim','var')
    clear Anim
end
Anim(nF) = struct('cdata',[],'colormap',[]);

% Time stamp
YLIM = get(gca,'ylim');
XLIM = get(gca,'xlim');

switch get(gca,'xdir')
    case 'normal'
        xloc = XLIM(1)+xloc*(XLIM(2)-XLIM(1));
    case 'reverse'
        xloc = XLIM(2)-xloc*(XLIM(2)-XLIM(1));
end
switch get(gca,'ydir')
    case 'normal'
        yloc = YLIM(1)+yloc*(YLIM(2)-YLIM(1));
    case 'reverse'
        yloc = yLIM(2)-yloc*(YLIM(2)-YLIM(1));
end

if RELATIVE_SB
    YLIM = get(hax,'ylim');
    XLIM = get(hax,'xlim');
    SB_X = XLIM(1)+SB_POS(1)*(XLIM(2)-XLIM(1)) + [0,PX_LENGTH];
    SB_Y = YLIM(1)+SB_POS(2)*(YLIM(2)-YLIM(1)) + [0,0];
end

for f=1:nF
    cla(hax);
    %hax = newplot(him)
    %make rgb data
    imbase = ind2rgb( gray2ind( mat2gray(origstack(:,:,f),clim),size(cmap,1)),cmap);
 
    %show rgb image
    set(him,'CData',imbase);
    
    
    axis(hax,'xy','image','off');
    hold(hax,'on');
    
    %Time Stamp
    if SHOW_TIMESTAMP
        str = sprintf('Time: %04.01f min',(CellData.Time(f)-CellData.Time(1))/60);
        text(xloc,yloc,str,'Color','w','VerticalAlignment','top','HorizontalAlignment','Left','FontSize',TIME_FONT_SIZE);
    end
    
    %Plot ScaleBar
    text(mean(SB_X),mean(SB_Y)+2*SB_WIDTH,sprintf('%d µm',SB_LENGTH),'parent',hax,'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
    plot(hax,SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);
    
    
    %% Plot KymoData
    if ~isempty(vel{f})
        %plot perimeter line
        s=linspace(ppF(f).breaks(1),ppF(f).breaks(end),500);
        xy = ppval(ppF(f),s);

        plot(hax,xy(1,:),xy(2,:),'-','Color','w','LineWidth',1.5);

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
        for j=1:numel(vel{f})
            quiver(hax,markersF{f}(1,j)',markersF{f}(2,j)',Vx(j),Vy(j),0,'Color',colors(j,:));
        end
    end
    
    

    %save frame to animation
    Anim(f) = getframe(hax);
end

%% Save mp4

[~,name] = fileparts(File);
[mov_file,mov_path] = uiputfile('*.mp4','Save movie?',fullfile(Dir,[name,'.mp4']));
if mov_file~=0
    writerObj = VideoWriter(fullfile(mov_path,mov_file),'MPEG-4');
    writerObj.FrameRate = 5;
    open(writerObj);
    writeVideo(writerObj,Anim);
    close(writerObj);
end



%% Play Movie
close(hfig);
implay(Anim);
putvar(Anim);


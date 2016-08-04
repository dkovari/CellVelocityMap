function PlotKymograph(CellData)
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
Lpts = CellData.kymo_Lpts;
vel = CellData.kymo_vel;
Time = CellData.Time;
Area = CellData.Area;
PxScale = CellData.PxScale;
L = CellData.kymo_L;


%% Load Colormaps
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
load(fullfile(pathstr,'KymoFunctions','SymetricVelcolormap.mat'),'mycmap')
load(fullfile(pathstr,'KymoFunctions','kymocmap.mat'),'kymocmap');
cmap = mycmap;

%% Fig Init
hfig = figure();
set(gcf,'Units','inches');
set(gcf,'Position',[1,1,5,8]);  %set figure size

%% Plot Kymograph
hAx = subplot(7,1,3:5);
velLim = [-0.15,0.15]; %color limits µm/sec

vh = 0;
for f=1:size(vel,1)
    vh = max(numel(Lpts{f}),vh);
end
vc = fix(vh/2)+1;
Vimg = NaN(vh,size(vel,1));
for f=1:size(vel,1)
    if isempty(vel{f})
        continue;
    end
    nl = numel(Lpts{f});
    nl2 = fix(nl/2);
    Vimg(vc:vc+nl2,f) = vel{f}(1:nl2+1);
    Vimg(vc-nl2:vc-1,f)=vel{f}(nl2+2:end);
end
Vimg(isnan(Vimg)) = 0;
cla(hAx,'reset');
%Vimg = radial_blur(Vimg,1,'method','gaussian');
hIm = imagesc(Vimg);
%set(hIm,'AlphaData',~isnan(Vimg));

%setup image width/height
l_low = Inf;
l_high = -Inf;
for f=1:size(vel,1)
    if isempty(vel{f})
        continue;
    end
    %Lpts{f} = Lpts{f};
    l_low = min(min(Lpts{f}),l_low);
    l_high = max(max(Lpts{f}),l_high);
end
set(hIm,'YData',[l_low,l_high]*PxScale);
axis xy;
set(hIm,'XData',[Time(1)-Time(1),Time(end-1)-Time(1)]/60);
axis auto

colormap(cmap);
set(gca,'CLim',velLim);
hold on;
% for f=1:size(vel,1)
%     plot([Time(f)/60-.5,f+.5],[-L(f)/2,-L(f)/2]*PxScale,'-k');
%     plot([f-.5,f+.5],[L(f)/2,L(f)/2]*PxScale,'-k');
% end

plot((Time(1:end-1)-Time(1))/60,L/2*PxScale,'-k','LineWidth',2);
plot((Time(1:end-1)-Time(1))/60,-L/2*PxScale,'-k','LineWidth',2);
ylabel('Perimeter [µm]','FontSize',12);
xlabel('Time [min]','FontSize',12);
hcb = colorbar;
hcbtitle = get(hcb,'ylabel');
set(hcbtitle,'string','Normal Velocity [µm/s]','FontSize',12);
%colorbar ticks
nCTicks = 5;
CBticks = linspace(velLim(1) ,velLim(2) ,nCTicks);
set(hcb,'YLimMode','manual');
set(hcb,'YLim',velLim);
set(hcb,'YTickMode','manual');
set(hcb,'YTick',linspace(velLim(1) ,velLim(2),nCTicks));
set(hcb,'YTickLabel',num2str(CBticks','%0.2f'));

xlim([Time(1)-Time(1),Time(end-1)-Time(1)]/60);

set(hAx,...
    'Box', 'off',...
    'FontSize', 12,...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Visible','on',...
    'TickDir','out',...
    'LineWidth',1.5);

%% Plot area time curve
hAx2 = subplot(7,1,1:2);
cla(hAx2,'reset');
plot(hAx2,(Time-Time(1))/60,Area,'-r','LineWidth',1.5);
set(hAx2,...
    'Box', 'off',...
    'FontSize', 12,...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Visible','on',...
    'TickDir','out',...
    'LineWidth',1.5);
xlim([Time(1)-Time(1),Time(end-1)-Time(1)]/60);
ylabel('Area [µm^2]','FontSize',12);
%xlabel('Time [min]','FontSize',12);

%make axes same width as kymograph by showing the colorbar
colorbar;

%% Plot Avearged Velocity
hAx3 = subplot(7,1,6:7);
cla(hAx3,'reset');
%zero line
hold on;plot(hAx3,[0,Time(end-1)-Time(1)],[0,0],':k','LineWidth',0.5);

AvgVel = cellfun(@nanmean,vel);
plot(hAx3,(Time(1:end-1)-Time(1))/60,AvgVel,'-k','LineWidth',1.5);
set(hAx3,...
    'Box', 'off',...
    'FontSize', 12,...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Visible','on',...
    'TickDir','out',...
    'LineWidth',1.5);
xlim([Time(1)-Time(1),Time(end-1)-Time(1)]/60);
ylabel('Avg. Velocity [µm/s]','FontSize',12);
xlabel('Time [min]','FontSize',12);
%make axes same width as kymograph by showing the colorbar
colorbar;

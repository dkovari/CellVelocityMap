function CalculateKymograph()

%add sub-functions to path
fpath = mfilename('fullpath');
pathstr = fileparts(fpath);
addpath(genpath(fullfile(pathstr,'KymoFunctions')));

%% Load Data
persistent last_dir;
%select file
[File,Dir] = uigetfile(fullfile(last_dir,'*.mat'),'Select cell data file');
if File==0
    return
end
if ~isempty(Dir)
    last_dir = Dir;
end

CellData = load(fullfile(Dir,File));

%% Validate File
if ~isfield(CellData,'threshstack')
    error('Data file does not contain threshstack variable');
end
if ~isfield(CellData,'Time')
    error('Data file does not contain Time variable');
end
if ~isfield(CellData,'origstack')
    error('Data file does not contain origstack variable');
end

if any(size(CellData.threshstack)~=size(CellData.origstack))
    error('origstack and threshstack are not the same size');
end
if numel(CellData.Time)~=size(CellData.threshstack,3)
    error('Number of time points doesnt match number of frames in images');
end

%% Check for other variable add if needed
if ~isfield(CellData,'PxScale')
    CellData.PxScale = 0.157825;
    save(fullfile(Dir,File),'-struct','CellData','PxScale','-append');
end

if ~isfield(CellData,'kymo_dL')
    CellData.kymo_dL = 3;
    save(fullfile(Dir,File),'-struct','CellData','kymo_dL','-append');
end

if ~isfield(CellData,'Area')
    CellData.Area = zeros(numel(CellData.Time),1);
    for f=1:numel(CellData.Time)
        CellData.Area(f) = nansum(nansum(CellData.threshstack(:,:,f)))*CellData.PxScale^2;
    end
    save(fullfile(Dir,File),'-struct','CellData','Area','-append');
end

%% Smooth Threshold Function
stack = largestBWstackregion(CellData.threshstack);
[H,W,nF] = size(stack);

if H*W>360*360
    RESIZE = true;
    data = radial_blur_stack(imresizestack(stack,0.5),5,'method','gaussian')-0.5;
%     figure(); stackfig(threshstack);colormap gray;
%     hold on; plot(line_pos(:,1)/2,line_pos(:,2)/2);
%     drawnow;
else
    RESIZE = false;
    data = radial_blur_stack(stack,5,'method','gaussian')-0.5;
end
if ~RESIZE
    perim = bwperimstack(data>=0);
else
    perim = bwperimstack(imresizestack(data,[H,W])>=0);
end
perimOL = imoverlaystack(CellData.origstack,perim,'Color',[1,1,0]);


%% Zero Line
if ~isfield(CellData,'kymo_zeroline')
    answer = questdlg('Loaded data does not specify zero line. Do you want to select one?','Zero Line');
    if ~strcmpi('yes',answer);
        return;
    end
    GET_ZL = true;
else
    answer = questdlg(...
        sprintf('Found Zero Line: [%d,%d]->[%d,%d]. Do you want to select a new one?',...
            CellData.kymo_zeroline(1,1),...
            CellData.kymo_zeroline(1,2),...
            CellData.kymo_zeroline(2,1),...
            CellData.kymo_zeroline(2,2)),...
        'Zero Line');
    if strcmpi('cancel',answer)
        return
    end
    if strcmpi('yes',answer)
        GET_ZL = true;
    else
        GET_ZL = false;
    end
end
if GET_ZL
    [xy0,xy1] = stack_line(perimOL);
    CellData.kymo_zeroline = [xy0;xy1];
    answer = questdlg('Save Zero Line to file?');
    if strcmpi('yes',answer)
        save(fullfile(Dir,File),'-struct','CellData','kymo_zeroline','-append');
    end
end

%% Start and End Time
if ~isfield(CellData,'kymo_start')
    CellData.kymo_start = 1;
end
if ~isfield(CellData,'kymo_end')
    CellData.kymo_end = nF;
end

hFig = stackfig(perimOL);
hold on;
plot(CellData.kymo_zeroline(:,1),CellData.kymo_zeroline(:,2),'-r');
title('Find Start/End Frames');

hMsg = msgbox('Find start and end frames, then hit ok to enter them in next window');

while ishandle(hMsg)
    pause(0.01);
    drawnow;
end
close(hFig);
drawnow;

prompt = {'Start Frame','End Frame'};
def = {num2str(CellData.kymo_start),num2str(CellData.kymo_end)};
while true
    answer = inputdlg(prompt,'Start/End Frames',1,def);
    if ~isempty(answer)
        k_s = str2double(answer{1});
        k_e = str2double(answer{2});
        if isnan(k_s)||isnan(k_e)
            continue;
        else
            CellData.kymo_start = max(1,k_s);
            CellData.kymo_end = min(nF,k_e);
            save(fullfile(Dir,File),'-struct','CellData','kymo_start','kymo_end','-append');
            break
        end
    else
        return;
    end
end


%% Calculate Kymograph
dT = diff(CellData.Time);
dL = CellData.kymo_dL;
PxScale = CellData.PxScale;
vel = cell(nF-1,1);
Lpts = cell(nF-1,1);
L = NaN(nF-1,1);
markersF = cell(nF-1,1);
ppF(nF-1) = struct('form','pp','breaks',[],'coefs',[],'pieces',0,'order',0,'dim',0);
ssF = cell(nF-1,1);

line_pos = CellData.kymo_zeroline;

%% loop over frames and compute velocity data
for f=CellData.kymo_start:CellData.kymo_end-1
    fprintf('Frame: %d/%d\n',f,CellData.kymo_end-1);
    %figure(99); imagesc(data(:,:,f))
    if all(all(data(:,:,f)==-0.5))||...
            all(all(isnan(data(:,:,f))))||...
            all(all(data(:,:,f)==1))||...
            all(all(data(:,:,f)==0))||...
            CellData.Area(f)==0||...
            isnan(CellData.Area(f))||...
            all(all(data(:,:,f+1)==-0.5))||...
            all(all(isnan(data(:,:,f+1))))||...
            all(all(data(:,:,f+1)==1))||...
            all(all(data(:,:,f+1)==0))||...
            CellData.Area(f+1)==0||...
            isnan(CellData.Area(f+1))
        fprintf('\tSkipping Frame %3\n',f);
        continue;
    end
    %% Use LS toolbox to interpolate intermediate steps
    LSstack = LSnormalevolve(data(:,:,f),data(:,:,f+1),...
                                'linearize','none',...
                                'nSteps',50,...
                                'speedFn','tanh06',...
                                'OutputStack',true,...
                                'timeStep',.1,...
                                'residual',10^-5);
    LSstack = cat(3,data(:,:,f),LSstack,data(:,:,f+1));
    
    if RESIZE
        LSstack = imresizestack(LSstack,[H,W]);
    end
    
    %% calculate initial markers
    %Convert first frame to explicit function
    pp=implicit2explicit(LSstack(:,:,1),'LargestOnly',true);
    %calc length
    [cumL,SL]= pplength_linear(pp);
    ppL = csapi(SL,cumL);
    ppS = csapi(cumL,SL);
    L(f) = cumL(end);
    %nMarkers = 2*fix(L(f)/dL/2)+1;

    %find the "zero" marker 
    So=intersectSpline2d(pp,[pp.breaks(1),pp.breaks(end)],line_pos(:,1),line_pos(:,2));
    if isempty(So)
        error('So empty')
    end
    if numel(So)>1, %if there are more than one intersections, get the one closest to the end
        xy = ppval(pp,So);
        dist = sqrt( (line_pos(2,1)-xy(1,:)).^2 + (line_pos(2,2)-xy(2,:)).^2);
        [~,mind] = min(dist);
        So = So(mind);
    end
    Lo = ppval(ppL,So);
    nM2 = fix(L(f)/dL/2);
    Lpts{f} = [ 0,(1:nM2)*dL,dL*(-nM2:1)];
    Lp = mod(Lo+Lpts{f},L(f));
    markersSS = ppval(ppS,Lp)';
    markersSS = pp_periodic(markersSS,pp);


    initMarkers.PP = pp;
    initMarkers.SS = markersSS;

    %% Starting with the precalculated markers, propegate markers to next
    %frame
    [markers,pp,ss] = propagateLSMarkers(LSstack,'initMarkers',initMarkers);
    markersF{f} =markers(:,:,1);
    ppF(f) = pp(1);
    ssF{f} = ss(:,1);

    %calculate distance moved along each marker
    m2 = markers(:,:,2:end);
    m1 = markers(:,:,1:end-1);
    dist = sqrt( (m2(1,:,:)-m1(1,:,:)).^2 + (m2(2,:,:)-m1(2,:,:)).^2 );
    dist = squeeze(dist);
    cumdist = cumsum(dist,2);

    vel{f} = cumdist(:,end)/dT(f)*PxScale;
    %vel(:,f-start_frame+1) = cumdist(:,end)/dT;

    %get sign of velocity
    [Gx,Gy] = gradient(LSstack(:,:,1));
    Vx = interp2(Gx,markers(1,:,1),markers(2,:,1));
    Vy = interp2(Gy,markers(1,:,1),markers(2,:,1));

    Dx = markers(1,:,end)-markers(1,:,1);
    Dy = markers(2,:,end)-markers(2,:,1);

    d = dot([Dx;Dy],[-Vx;-Vy]);
    vel{f} = vel{f}.*sign(d)';
end

if RESIZE
    data = imresizestack(data,[H,W]);
end

%% Replace Errors with NaN
if ~isfield(CellData,'kymo_velLimit')
    CellData.kymo_velLimit = 0.3;
end
while true
    answer = inputdlg({'Velocity Limit'},'Velocity Limit',1,{num2str(CellData.kymo_velLimit)});
    if isempty(answer)
        break
    end
    v = str2double(answer{1});
    if ~isnan(v)&&v>0
        CellData.kymo_velLimit = v;
        break;
    end
end
save(fullfile(Dir,File),'-struct','CellData','kymo_velLimit','-append');
for f=1:numel(vel)
    vel{f}(abs(vel{f})>CellData.kymo_velLimit) = NaN;
end

%% Save KymoData
CellData.kymo_Lpts = Lpts;
CellData.kymo_vel = vel;
CellData.kymo_markersF = markersF;
CellData.kymo_data = data;
CellData.kymo_ppF = ppF;
CellData.kymo_ssF = ssF;
CellData.kymo_L = L;

save(fullfile(Dir,File),'-struct',...
    'CellData',...
    'kymo_Lpts',...
    'kymo_vel',...
    'kymo_markersF',...
    'kymo_data',...
    'kymo_ppF',...
    'kymo_ssF',...
    'kymo_L',...
    '-append');

%% Plot Kymograph
PlotKymograph(CellData);
AnimateKymograph(CellData,fullfile(Dir,File));

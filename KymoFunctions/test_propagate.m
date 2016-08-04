%% Test PrrpagetLSMarkers

%% Generate fake data

% data = zeros(256,256,25);
% for f=1:25
%     data(:,:,f) = double(BWcircle(10+f*2,[128,128],[256,256]));
% end
% data = radial_blur_stack(data,10,'method','gaussian')-0.5;

% x0 = 128; y0= 128;
% %line which defines zero-th markers
% Xo = [x0,256]';
% Yo = [y0,y0]';
% 
% d0 = double(BWcircle(25,[128,128],[256,256]));
% d1 = zeros(256);
% d1(32:224,32:224) = 1;

%% use real data, must preload a cell movie stack (use track reproc)
load('..\..\KymoData_old\2014-05-18-belb-s1-t7-stackdata.mat')
threshstack = largestBWstackregion(threshstack);
d0 = threshstack(:,:,30);
d1 = threshstack(:,:,37);

%Get initial centroid position
stats = regionprops(d0,'Centroid');
x0 = stats.Centroid(1);
y0 = stats.Centroid(2);

[H,W] = size(d0);

%line which defines zero
Xo = [x0,W]';
Yo = [y0,y0]';


%% blur the data to smooth edge and create a continious implicit function
d0 = radial_blur(d0,10,'method','gaussian')-0.5;
d1 = radial_blur(d1,10,'method','gaussian')-0.5;

%% reinitialize the implicit function so that the gradient is always 1 near the 0th levelset
%see O&F for reason why
%this step is not entirely necessary if the bluring function does a good
%job creating a constant gradient near edge
d0 = reinit(d0,'NumCells',30);
d1 = reinit(d1,'NumCells',30);

%% Use LS toolbox to interpolate intermediate steps
[data,t,r] = LSnormalevolve(d0,d1,...
                            'linearize','none',...
                            'nSteps',300,...
                            'speedFn','tanh08',...
                            'OutputStack',true,...
                            'timeStep',.1,...
                            'residual',10^-4);
data = cat(3,d0,data,d1);


%% setup initial markers
%Convert first frame to explicit function
pp=implicit2explicit(data(:,:,1),'LargestOnly',true);
%calc length
[cumL,SL]= pplength_linear(pp);
ppL = csapi(SL,cumL);
ppS = csapi(cumL,SL);
L = cumL(end);
dL = 5;
nMarkers = 2*fix(L/dL/2)+1;

%find the "zero" marker 
So=intersectSpline2d(pp,[pp.breaks(1),pp.breaks(end)],Xo,Yo);
if isempty(So)
    error('So empty')
end
Lo = ppval(ppL,So(1));
nM2 = fix(L/dL/2);
Lpts = [ 0,(1:nM2)*dL,dL*(-nM2:1)];
Lp = mod(Lo+Lpts,L);
markersSS = ppval(ppS,Lp)';
markersSS = pp_periodic(markersSS,pp);


initMarkers.PP = pp;
initMarkers.SS = markersSS;

%% propagate markers
[markers,pp,ss] = propagateLSMarkers(data,'initMarkers',initMarkers,'DisplayDebug',false);

%% generate plot
figure();
hold on;
imagesc(origstack(:,:,30));
colormap gray;
axis image;
axis xy;
axis off;
for f = 1:size(markers,3)-1:size(markers,3)
    %plot edge line
    xy = ppval(pp(f),linspace(pp(f).breaks(1),pp(f).breaks(end),1000));
    plot(xy(1,:),xy(2,:),'-y');
end

%plot markers with path
for m=1:size(markers,2)
    plot(reshape(markers(1,m,:),[],1),reshape(markers(2,m,:),[],1),'.-r');
end

% set limits
xlim([120,240]);
ylim([120,240]);

%plot scale bar
PxScale = 0.157825; %um/px

SBL = 125;
SBR = 5/PxScale+SBL;
SBY = 125;
plot([SBL,SBR],[SBY,SBY],'-w','LineWidth',3);

text((SBL+SBR)/2,SBY+1,'5 µm',...
    'Color','w',...
    'FontSize',18,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom');

%fig2pdf('LSDemo_interp.pdf')

%% Propagate markers without levelset interpolation
[m1,pp1,ss1] = propagateLSMarkers(cat(3,d0,d1),'initMarkers',initMarkers);
figure();
hold on;
imagesc(origstack(:,:,30));
colormap gray;
axis image;
axis xy;
axis off;
for f = 1:size(m1,3)-1:size(m1,3)
    %plot edge line
    xy = ppval(pp1(f),linspace(pp1(f).breaks(1),pp1(f).breaks(end),1000));
    plot(xy(1,:),xy(2,:),'-y');
end

%plot markers with path
for m=1:size(m1,2)
    plot(reshape(m1(1,m,:),[],1),reshape(m1(2,m,:),[],1),'.-r');
end

% set limits
xlim([120,240]);
ylim([120,240]);

%plot scale bar
PxScale = 0.157825; %um/px

SBL = 125;
SBR = 5/PxScale+SBL;
SBY = 125;
plot([SBL,SBR],[SBY,SBY],'-w','LineWidth',3);

text((SBL+SBR)/2,SBY+1,'5 µm',...
    'Color','w',...
    'FontSize',18,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom');

%fig2pdf('LSDemo_nointerp.pdf')
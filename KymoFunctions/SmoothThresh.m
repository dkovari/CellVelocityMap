%% Reprecess the cell image to create a smooth threshold

STD = Inf;
ILOW = -0.003;
IHIGH = 0.002;

th2 = stackthresh(origstack,STD,[ILOW,IHIGH]);
th2 = th2&threshstack;
th2 = imfillstack(th2,'holes');
th2=bwareaopenstack(th2,16);
th2 = radial_blur_stack(th2,10,'method','gaussian')>0.45;
th2 = radial_blur_stack(th2,10,'method','gaussian')>0.55;
th2 = imfillstack(th2,'holes');
%% View Stack
OLs = zeros(size(th2));
for f=1:size(th2,3)
    OL = bwperim(th2(:,:,f));
    im = origstack(:,:,f);
    im(OL) = NaN;
    OLs(:,:,f) = im;
end
stackfig(OLs,'figure',figure(2),'clim','average'); colormap gray;
% userdata.Time = Time;
% stackfig(OLs,'figure',figure(2),'clim','average','frameupdate_fn',@show_text_fn,'userdata',userdata); colormap gray;

Area = PxScale^2*squeeze(sum(sum(th2,1),2));

figure(4); plot(Time,Area,'-r');

if exist('FILE','var')
    save(fullfile(DIR,FILE),'th2','-append');
end

clear OL OLs im
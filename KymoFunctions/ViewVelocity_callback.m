function ViewVelocity_callback(hObject)

handles = guidata(hObject);
f = handles.curFrame;
if ~isfield(handles,'hTxt_time')||~ishghandle(handles.hTxt_time)
    xlim = get(handles.hAx_Img,'xlim');
    ylim = get(handles.hAx_Img,'ylim');
    xpos = xlim(1)+0.03*(xlim(2)-xlim(1));
    ypos = ylim(2)-0.03*(ylim(2)-ylim(1));
    
    handles.hTxt_time = text(xpos,ypos,'Time:',...
        'parent',handles.hAx_Img,...
        'Color','w',...
        'VerticalAlignment','Cap',...
        'HorizontalAlignment','Left',...
        'FontSize',14);
end

set(handles.hTxt_time,'String',sprintf('Time: %06.1f',handles.userdata.Time(f)));


if f>handles.nFrames-1
    
    if ~isfield(handles,'hVecs')
        handles.hVecs = [];
    end
    delete(handles.hVecs(ishghandle(handles.hVecs)));
    
    guidata(hObject,handles);
    if ~isfield(handles,'hOL')||(~isempty(handles.hOL)&&~ishghandle(handles.hOL))
        handles.hOL = [];
    end
    delete(handles.hOL);
    guidata(hObject,handles);
    return;
end




%plot outline
if ~isempty(handles.userdata.vel{f})
    % generate line data
    s=linspace(handles.userdata.ppF(f).breaks(1),handles.userdata.ppF(f).breaks(end),500);
    xy = ppval(handles.userdata.ppF(f),s);
    
    if ~isfield(handles,'hOL')||(~isempty(handles.hOL)&&~ishghandle(handles.hOL))
        handles.hOL = line('Parent',handles.hAx_Img,'XData',xy(1,:),'YData',xy(2,:),'Color','y');
    else
        set(handles.hOL,'XData',xy(1,:),'YData',xy(2,:));
    end

    %vector data
    [Gx,Gy] = gradient(handles.userdata.data(:,:,f));
    Vxg = interp2(Gx,handles.userdata.markersF{f}(1,1)',handles.userdata.markersF{f}(2,1)');
    Vyg = interp2(Gy,handles.userdata.markersF{f}(1,1)',handles.userdata.markersF{f}(2,1)');

    %calp perp line
    dxy = fnval(fnder(handles.userdata.ppF(f),1),handles.userdata.ssF{f})';
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
    Vx = Vx.*handles.userdata.vel{f}*handles.userdata.VSCALE;
    Vy = Vy.*handles.userdata.vel{f}*handles.userdata.VSCALE;

    %if needed delete old vectors
    if ~isfield(handles,'hVecs')
        handles.hVecs = [];
    end

    delete(handles.hVecs(ishghandle(handles.hVecs)));
    %% plot vectors
    %color data
    colors = num2climcolor(handles.userdata.mycmap,handles.userdata.climVEL,handles.userdata.vel{f});
    nV = numel(handles.userdata.vel{f});

    %plot last
    X = [handles.userdata.markersF{f}(1,nV)', Vx(nV)+handles.userdata.markersF{f}(1,nV)'];
    Y = [handles.userdata.markersF{f}(2,nV)', Vy(nV)+handles.userdata.markersF{f}(2,nV)'];
    handles.hVecs(nV) = line('Parent',handles.hAx_Img,'XData',X,'YData',Y,'Color',colors(nV,:));
    for n=1:nV-1
        X = [handles.userdata.markersF{f}(1,n)', Vx(n)+handles.userdata.markersF{f}(1,n)'];
        Y = [handles.userdata.markersF{f}(2,n)', Vy(n)+handles.userdata.markersF{f}(2,n)'];
        handles.hVecs(n) = line('Parent',handles.hAx_Img,'XData',X,'YData',Y,'Color',colors(n,:));
    end
else
    %delete line
    if ~isfield(handles,'hOL')||(~isempty(handles.hOL)&&~ishghandle(handles.hOL))
        handles.hOL = [];
    end
    delete(handles.hOL);
    %delete vectors
    if ~isfield(handles,'hVecs')
        handles.hVecs = [];
    end
    delete(handles.hVecs(ishghandle(handles.hVecs)));
end

guidata(hObject,handles);
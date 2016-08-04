function RGB = num2climcolor(cmap,clim,val)

if ischar(cmap)
    cmap = eval([cmap,'(256)']);
    if isempty(cmap)
        error('cmap must be a valid colormap');
    end
end

if clim(1)>=clim(2)
    error('clim must be an increasing range')
end

%columnize val
val = val(:);
val(val<clim(1)) = clim(1);
val(val>clim(2)) = clim(2);
val = floor( (val-clim(1))/(clim(2)-clim(1))*(size(cmap,1)-1))+1;
val(isnan(val)) = 1;
RGB = cmap(val,:);

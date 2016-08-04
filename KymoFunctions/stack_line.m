function [xy0,xy1] = stack_line(stack)
%display stackfig and allow user to pick a line.
%returns coorinates of line after user closes figure

hf = stackfig(stack,'Colormap','gray');
hax = get(hf,'CurrentAxes');


[H,W,~] = size(stack);

xy0 = [W,H]/2;
xy1 = [W-2,H-2];

hl = imline(hax,[xy0;xy1]);

title(hax,'Choose line, close figure to continue')

%%inline function
    function return_line(~,~)
        pos = hl.getPosition();
        xy0 = pos(1,:);
        xy1 = pos(2,:);
        delete(hf);
    end

set(hf,'CloseRequestFcn',@return_line);

waitfor(hf);

end

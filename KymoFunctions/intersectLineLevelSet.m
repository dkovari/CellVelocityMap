function [X,Y,C,pp,T] = intersectLineLevelSet(data,xy0,dxy,varargin)
%find the intersection of a line with a levelset of the surface defined by
%gridded data
% Inputs:
%   data: 2-D gridded data defining a z=f(x,y) surface
%   [x0,y0]: a point on the line
%   [dx,dy]: a vector along the line
% Parameters:
%   'LevelSet', value:  value of the levelset to find intersection with
%               default: 0
% Requirements:
%   coutourcs()
%   implicit2explicit()

p=inputParser;
p.KeepUnmatched=true;
p.CaseSensitive=false;
addParamValue(p,'LevelSet',0,@(x) isnumeric(x)&&isscalar(x));
parse(p,varargin{:});

LS = p.Results.LevelSet;

if size(xy0,1)~=size(dxy,1)
    error('size of xy0,dxy must be the same');
end

NLines = size(xy0,1);

% Xn=NaN(size(xy0,1),2);
% Yn=Xn;
% for n=1:size(xy0,1)
%     if dxy(n,1)==0    %vertical line spanning LS image
%         Xn(n,:)=[1,1]*xy0(n,1);
%         Yn(n,:)=[0,size(data,1)+1];
%     else
%         xn = [0,size(data,2)+1];
%         yn = (dxy(n,2)/dxy(n,1))*(xn-xy0(n,1))+xy0(n,2);
%         %set line limits to boundary of the LS image
%         [xn,yn]=intersections(xn,yn,[0,size(data,2)+1,size(data,2)+1,0,0],[0,0,size(data,1)+1,size(data,1)+1,0]);
%         Xn(n,:) = xn(:);
%         Yn(n,:) = yn(:);
%     end
% end

%find the countour of the levelset
[pp,~,C] = implicit2explicit(data,'LevelSet',LS,'LargestOnly',true,'ForcePeriodic',true);


debug = false;
if debug
    figure(2);
    clf;
    imagesc(data); colormap gray; axis image; axis xy; hold on;
    plot(C.X,C.Y,'.g');
end

%find breakpoints which cross line
m = dxy(:,2)./dxy(:,1);
Yfn = @(x) m*reshape(x,1,[])-repmat((m.*xy0(:,1)-xy0(:,2)),1,numel(x));

Ycx = Yfn(C.X);


CY = repmat(C.Y,NLines,1);
CX = repmat(C.X,NLines,1);

Vert = dxy(:,1)==0;
Vert = repmat(Vert,1,size(CX,2));
Xl = repmat(xy0(:,1),1,size(CX,2));

B = (CY(:,1:end-1)<=Ycx(:,1:end-1)&CY(:,2:end)>=Ycx(:,2:end))|...
    (CY(:,1:end-1)>=Ycx(:,1:end-1)&CY(:,2:end)<=Ycx(:,2:end))|...
    Vert(:,1:end-1)&(...
    CX(:,1:end-1)<=Xl(:,1:end-1)&CX(:,2:end)>=Xl(:,2:end)|...
    CX(:,1:end-1)>=Xl(:,1:end-1)&CX(:,2:end)<=Xl(:,2:end));

X = cell(NLines,1);
Y = cell(NLines,1);
T = cell(NLines,1);

An = dxy(:,2);%(Yn(:,2)-Yn(:,1));
Bn = -dxy(:,1);%(Xn(:,1)-Xn(:,2));
Cn = xy0(:,1).*(-dxy(:,2))+xy0(:,2).*(dxy(:,1)); %Xn(:,1).*(-An)+Yn(:,1).*(-Bn);
[px,py] = ppsplitdim(pp);

for n=1:NLines
    Bidx = find(B(n,:));
    if ~isempty(Bidx)
        T{n} = pp.breaks(Bidx);
        T2 = pp.breaks(Bidx+1);
        for m=1:numel(Bidx)
            %create a new polynomial which solves the line-curve intersection
            poly = An(n)*px.coefs(Bidx(m),:)+Bn(n)*py.coefs(Bidx(m),:)+[zeros(1,px.order-1),Cn(n)];
            r = roots(poly);
            %filter r to get rid of imaginary and solutions outside of
            %break point limits
            tlim = T2(m)-T{n}(m);
            r(imag(r)~=0|r<-0.01|r>1.01*tlim) = [];
            r = r(1);   %this is a hack
            if numel(r)>1||isempty(r)
                r
                error('check r...dan needs to fix the polynomial solver');
            end
            if ~isempty(r)
                T{n}(m) = T{n}(m)+r(1);
            end
            if debug
                xy = ppval(pp,T{n});
                plot(xy(1,:),xy(2,:),'om');
                title(['n=',num2str(n)]);
            end
        end
        
        
%         T{n} = pp.breaks(Bidx);
%         T2 = pp.breaks(Bidx+1);
% %         xy = ppval(pp,T{n});
% %         plot(Xn(n,:),Yn(n,:),'-r');
% %         plot(xy(1,:),xy(2,:),'om');
% %         title(['n=',num2str(n)]);
% %         pause
%         for m=1:numel(T{n}) %loop over each break point and solve the roots of the polynomial
%             %T{n}(m) = fzero(@(t) distline(t,n),T{n}(m));
%             
%             
%             
%             
%             
%         end
    end
    xy = ppval(pp,T{n});
    X{n} = xy(1,:);
    Y{n} = xy(2,:);
    
%     if any(size(T{n})==0)
%         Bidx
%         Y{n}
%         size(T{n})
%         n
%         
%         figure(2);
%         clf;
%         imagesc(data); colormap gray; axis image; axis xy; hold on;
%         plot(C.X,C.Y,'.g');
%         plot(xy0(:,1),xy0(:,2),'.b')
%         quiver(xy0(n,1),xy0(n,2),dxy(n,1)*10,dxy(n,2)*10,0,'-r');
%         %error('no S');
%     end
end

if debug
    for n=1:size(X,1)
        plot(X{n},Y{n},'xy');
    end
end

%% Define subfunction
% function d = distline(t,n)%pp,x1,x2,y1,y2)
% t=mod(t,max(pp.breaks));    %make t periodic
% fxy = ppval(pp,t);
% d = ( (Yn(n,2)-Yn(n,1))*fxy(1,:) -...
%     (Xn(n,2)-Xn(n,1))*fxy(2,:)...
%     + Xn(n,2)*Yn(n,1) - Yn(n,2)*Xn(n,1) )...
%     /sqrt( (Yn(n,2)-Yn(n,1))^2 + (Xn(n,2)-Xn(n,1))^2 );
% end




end




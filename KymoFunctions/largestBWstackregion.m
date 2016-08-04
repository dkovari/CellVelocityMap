function bwstack2 = largestBWstackregion(bwstack,conn)

if ndims(bwstack)>3
    error('bwstack must be 3d, [H,W,nFrames]')
end

if nargin<2
    conn = 8;
end

bwstack2 = bwstack;

for f = 1:size(bwstack,3)
    bwstack2(:,:,f) = largestBWregion(bwstack(:,:,f),conn);
end


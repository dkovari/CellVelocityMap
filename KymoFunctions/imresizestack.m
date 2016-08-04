function stack2 = imresizestack(stack,varargin)

im1 = imresize(stack(:,:,1),varargin{:});
stack2(size(im1,1),size(im1,2),size(stack,3)) = 0;
stack2(:,:,1) = im1;
for f=2:size(stack,3)
    stack2(:,:,f) = imresize(stack(:,:,f),varargin{:});
end
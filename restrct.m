function [fc] = restrct(ff,fc,nnx,nny,ncx,ncy)

f = ff;
% coarsening in both x and y
for ic=2:ncx-1
    for jc=2:ncy-1
        i = 2*ic-1;
        j = 2*jc-1;
        fc(ic,jc) = (4*f(i,j) + 2*(f(i+1,j) + f(i,j+1)+ ...
            f(i-1,j) + f(i,j-1)) + (f(i-1,j-1) + f(i+1,j-1) + f(i-1,j+1) + f(i+1,j+1)))/16;
    end
end
%
% Set residual on boundaries
%
ix = 1;
jy = 1;

for jc=1:ncy-1:ncy
    j = jc + jy*(jc-1);
    jm1= max(j-1,2);
    jp1 = min(j+1,nny-1);
    
    % Corners of the grid
    for ic=1:ncx-1:ncx
        i = ic + ix*(ic-1);
        im1 = max(i-1,2);
        ip1 = min(i+1,nnx-1);
        fc(ic,jc) = (4*f(i,j) + 2*(f(ip1,j) + f(i,jp1)+ ...
            f(im1,j) + f(i,jm1)) + (f(im1,jm1) + f(ip1,jm1) + f(im1,jp1) + f(ip1,jp1)))/16;
    end
    % y = yc,yd interior edges
    for ic=2:ncx-1
        i=2*ic-1;
        fc(ic,jc) = (4*f(i,j) + 2*(f(i+1,j) + f(i,jp1)+ ...
            f(i-1,j) + f(i,jm1)) + (f(i-1,jm1) + f(i+1,jm1) + f(i-1,jp1) + f(i+1,jp1)))/16;
    end
end

%     set x=xa,xb interior edges
for ic=1:ncx-1:ncx
    i = ic + ix*(ic-1);
    im1 = max(i-1,2);
	ip1 = min(i+1,nnx-1);
    for jc=2:ncy-1
        j = 2*jc-1;
        fc(ic,jc) = (4*f(i,j) + 2*(f(ip1,j) + f(i,j+1)+ ...
            f(im1,j) + f(i,j-1)) + (f(im1,j-1) + f(ip1,j-1) + f(im1,j+1) + f(ip1,j+1)))/16;
    end
end

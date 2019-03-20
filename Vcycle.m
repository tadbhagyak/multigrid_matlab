function [Tcorr] = Vcycle(T_fg,nxk,nyk,kcycle)
%% Set up V cycle of multigrid
% rs  o        o (Ts + Ih*deltah)
%      .      .
%       o    o
%        .  .
% deltah  x
%
rs_k = rs(2:nnx+1,1:nny);

% Restriction operator
ff = rs_k;
cur_grid = kcycle+1;

for l = cur_grid:kcycle
    
    nnx = nxk(l);
    nny = nyk(l);
    
    % Number of points on coarsest grid
    ncx = nxk(l-1);
    ncy = nyk(l-1);
    resc = zeros(ncx,ncy);
    
    % restrict the residual
    resc = restrct(ff,resc,nnx,nny,ncx,ncy); % (ncx,ncy)
    res1 = resc;
    
    % Regrid and find points on the new grid
    [xc,yc,etac,detadyc,der1c,L] = gridCluster(ncx,ncy);
    
    % Redefine material matrices on new grid
    [Ac,muc,psic] = materialMatrix(xc,L,ncx,ncy); %(ncx+2,ncy+2)
    
    % Coefficient stencil on new grid
    phic = surfTop(xc,L,ncx);
    [nwc,nc,nec,wc,cec,ec,sec,sc,swc] = coeff_stencil(Ac,muc,phic,der1c,detadyc,xc,yc,ncx,ncy);
    
    % Smooth the residual / Solve for delta-H
    niter = 10;
    if (l == 2) niter = 10000; end
    [deltac] = smoothRes(resc,ncx,ncy,nwc,nc,nec,ec,cec,wc,swc,sc,sec,niter);
    ff = zeros(ncx,ncy);
    
    % Calculate residual for deltaH - recursion here
    for i=2:ncx+1
        for j=2:ncy
            ic = i-1;
            ip = i+1;
            im = i-1;
            if (i==2) im = ncx; end
            if (i==ncx+1) ip = 3; end
            ff(ic,j) = resc(ic,j) -(nc(i,j)*deltac(i,j+1) + sc(i,j)*deltac(i,j-1) + cec(i,j)*deltac(i,j)...
                + wc(i,j)*deltac(im,j) + ec(i,j)*deltac(ip,j) + nec(i,j)*deltac(ip,j+1)...
                + nwc(i,j)*deltac (im,j+1) + sec(i,j)*deltac(ip,j-1) + swc(i,j)*deltac(im,j-1));
        end
    end
end
delta_cg = deltac(2:ncx+1,1:ncy);
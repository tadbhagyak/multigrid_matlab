function [Tcorr,delta,r_cg] = Vcycle(res_fg,nxk,nyk,cur_grid,niter)
%% Set up V cycle of multigrid
% rs  o        o (Ts + Ih*deltah)
%      .      .
%       o    o
%        .  .
% deltah  x
% 

% Restriction operator
ff = res_fg;
% arrays to store solution and residual at each level
delta_1 = [];
delta_2 = [];
delta_3 = [];
delta_4 = [];
r_cg1 = [];
r_cg2 = [];
r_cg3 = [];
r_cg4 = [];

for l = cur_grid:-1:2
    nnx = nxk(l);
    nny = nyk(l);
    
    % Number of points on coarsest grid
    ncx = nxk(l-1);
    ncy = nyk(l-1);
    resc = zeros(ncx,ncy);
    
    % restrict the residual
    resc = restrct(ff,resc,nnx,nny,ncx,ncy); %(ncx,ncy)
    res1 = resc;
    
    % Regrid and find points on the new grid
    [xc,yc,etac,detadyc,der1c,L] = gridCluster(ncx,ncy);
    
    % Redefine material matrices on new grid
    [Ac,muc,psic] = materialMatrix(xc,L,ncx,ncy); %(ncx+2,ncy+2)
    
    % Coefficient stencil on new grid
    phic = surfTop(xc,L,ncx);
    [nwc,nc,nec,wc,cec,ec,sec,sc,swc] = coeff_stencil(Ac,muc,phic,der1c,detadyc,xc,yc,ncx,ncy);
    
    % Smooth the residual / Solve for delta-H
    niter = 1;
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
    delta_cg = deltac(2:ncx+1,1:ncy);
    
    d_temp = reshape(delta_cg,ncx*ncy,1);
    f_temp = reshape(ff,ncx*ncy,1);
    
    % Store the delta and residuals on grids
    if (l == 2)
        delta_1 = delta_cg;
        r_cg1  = ff;
    elseif (l==3)
        delta_2 = delta_cg;
        r_cg2  = ff;
    elseif (l==4)
        delta_3 = delta_cg;
        r_cg3 = ff;
    end
end

%% Coarse grid correction of the residual
% Prolongation operator
% Consider two meshes VH(fine) and Vh(coarse). The prolongation operator
% takes vector from coarse to fine grid using linear interpolation i.e.
%   I(h-H) : VH ---> Vh
niter = 1;
for m = 1:cur_grid
    cg = m+1;
    nnx = nxk(cg);
    nny = nyk(cg);
    
    % Number of points on coarsest grid
    ncx = nxk(m);
    ncy = nyk(m);
    deltaf = zeros(nnx,nny);
    
    if (m==1)
        deltaf = prolong(deltaf,delta_1,nnx,nny,ncx,ncy);
        delta_2(2:nnx+1,1:nny) = delta_2(2:nnx+1,1:nny) + deltaf;
        delta_2(1,1:nny) = delta_2(nnx,1:nny);
        delta_2(nnx+2,1:nny) = delta_2(3,1:nny);
        [delta_2,res_ps] = postSmooth(delta_2,nxk,nyk,cg,niter);
        
    elseif (m==2)
        deltaf = prolong(deltaf,delta_2,nnx,nny,ncx,ncy);
        delta_3(2:nnx+1,1:nny) = delta_3(2:nnx+1,1:nny) + deltaf;
        delta_3(1,1:nny) = delta_3(nnx,1:nny);
        delta_3(nnx+2,1:nny) = delta_3(3,1:nny);
        [delta_3,res_ps] = postSmooth(delta_3,nxk,nyk,cg,niter);
        
    elseif (m==3)
        deltaf = prolong(deltaf,delta_3,nnx,nny,ncx,ncy);
        delta_4(2:nnx+1,1:nny) = delta_4(2:nnx+1,1:nny) + deltaf;
        delta_4(1,1:nny) = delta_4(nnx,1:nny);
        delta_4(nnx+2,1:nny) = delta_4(3,1:nny);
        [delta_4,res_ps] = postSmooth(delta_4,nxk,nyk,cg,niter);
        
    end
end

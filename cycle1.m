function [Tf_ps,rs_ps] = cycle1(nxk,nyk,kcycle)
    nx = nxk(kcycle); %finest grid points in x
    ny = nyk(kcycle); %finest grid points in y
    T = zeros(nx+2,ny+2); Tn = T; res= T;
    niter = 1000;
    
    % Setup the coarsest grid
    psi = ones(nx+2,ny+2);
    [x,y,eta,detady,der1,L] = gridCluster(nx,ny);
    phi = surfTop(x,L,nx);
    % Stencil and Coefficients for the material matrix
    [A,mu,psi] = materialMatrix(x,L,nx,ny);
    [nw,n,ne,w,ce,e,se,s,sw] = coeff_stencil(A,mu,phi,der1,detady,x,y,nx,ny);
    
    % Call the relaxation solver
    [Tn,rn] = postSmooth(T,nx,ny,nw,n,ne,e,ce,w,sw,s,se,niter);
    Tc = Tn(2:nx+1,1:ny);
    
    % Interpolate the solution to the next finest grid
    Tf = interpSol(Tc,nxk,nyk,kcycle);
    
    % Post smooth on the same grid
    [x,y,eta,detady,der1,L] = gridCluster(nnx,nny);
    phi = surfTop(x,L,nnx);
    [A,mu,psi] = materialMatrix(x,L,nnx,nny);
    [nw,n,ne,w,ce,e,se,s,sw] = coeff_stencil(A,mu,phi,der1,detady,x,y,nnx,nny);
    T_temp = zeros(nnx+2,nny+2);
    T_temp(2:nnx+1,1:nny) = Tf;
    T_temp(1,:) = T_temp(nnx,:);
    T_temp(nnx+2,:) = T_temp(3,:);
    [Ts,rs] = postSmooth(T_temp,nnx,nny,nw,n,ne,e,ce,w,sw,s,se,1);
    
    %% Set up V cycle of multigrid
    % rs  o        o (Ts + Ih*deltah)
    %      .      .
    %       .    .
    %        .  .
    % deltah  x
    %
    rs_k = rs(2:nnx+1,1:nny);
    % Restriction operator
    ff = rs_k;
    cur_grid = kcycle+1;
    
    for l = cur_grid:kcycle+1
        
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
    
    %% Coarse grid correction of the residual
    % Prolongation operator
    % Consider two meshes VH(fine) and Vh(coarse). The prolongation operator
    % takes vector from coarse to fine grid using linear interpolation i.e.
    %   I(h-H) : VH ---> Vh
    niter = 1;
    for m = 1%:2
        nnx = nxk(m+1);
        nny = nyk(m+1);
        % Number of points on coarsest grid
        ncx = nxk(m);
        ncy = nyk(m);
        deltaf = zeros(nnx,nny);
        deltaf = prolong(deltaf,delta_cg,nnx,nny,ncx,ncy);
        
        Ts(2:nnx+1,1:nny) = Ts(2:nnx+1,1:nny) + deltaf;
        Ts(1,1:nny) = Ts(nnx,1:nny);
        Ts(nnx+2,1:nny) = Ts(3,1:nny);
        [Tf_ps,rf_ps] = postSmooth(Ts,nnx,nny,nw,n,ne,e,ce,w,sw,s,se,niter);
    end
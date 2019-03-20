function [Tf_ps,rs_ps] = cycle1(nxk,nyk,kcycle)
    nx = nxk(kcycle); %finest grid points in x
    ny = nyk(kcycle); %finest grid points in y
    T = zeros(nx+2,ny+2); Tn = T; res= T;
    niter = 1000;
       
    % Call the relaxation solver on the coarsest grid
    cur_grid = kcycle;
    [Tn,res] = postSmooth(T,nxk,nyk,cur_grid,niter);
    Tc = Tn(2:nx+1,1:ny);
    
    % Interpolate the solution to the next finest grid
    Tf = interpSol(Tc,nxk,nyk,kcycle);
    
    % Post smooth on the finer grid
    cur_grid = kcycle+1;
    nnx = nxk(cur_grid);
    nny = nyk(cur_grid);
    T_temp = zeros(nnx+2,nny+2);
    T_temp(2:nnx+1,1:nny) = Tf;
    T_temp(1,:) = T_temp(nnx,:);
    T_temp(nnx+2,:) = T_temp(3,:);
    [Ts_fg,res_fg] = postSmooth(T_temp,nxk,nyk,cur_grid,niter);
    
    %% Set up V cycle of multigrid
    % rs  o        o (Ts + Ih*deltah)
    %      .      .
    %       .    .
    %        .  .
    % deltah  x
    % Restriction operator
    cur_grid = kcycle+1;
    
    for l = cur_grid:2
        
        nnx = nxk(l);
        nny = nyk(l);
        
        % Number of points on coarsest grid
        ncx = nxk(l-1);
        ncy = nyk(l-1);
        resc = zeros(ncx,ncy);
        
        % restrict the residual
        resc = restrct(res_fg,resc,nnx,nny,ncx,ncy); % (ncx,ncy)
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
    end
    delta_cg = deltac(2:ncx+1,1:ncy);
    
    %% Coarse grid correction of the residual
    % Prolongation operator
    % Consider two meshes VH(fine) and Vh(coarse). The prolongation operator
    % takes vector from coarse to fine grid using linear interpolation i.e.
    %   I(h-H) : VH ---> Vh
    niter = 1;
    for m = 1%:2
        cur_grid = m+1;
        nnx = nxk(cur_grid);
        nny = nyk(cur_grid);
        % Number of points on coarsest grid
        ncx = nxk(m);
        ncy = nyk(m);
        deltaf = zeros(nnx,nny);
        deltaf = prolong(deltaf,delta_cg,nnx,nny,ncx,ncy);
        
        Ts_fg(2:nnx+1,1:nny) = Ts_fg(2:nnx+1,1:nny) + deltaf;
        Ts_fg(1,1:nny) = Ts_fg(nnx,1:nny);
        Ts_fg(nnx+2,1:nny) = Ts_fg(3,1:nny);
        [Tf_ps,rs_ps] = postSmooth(Ts_fg,nxk,nyk,cur_grid,niter);
        %Tf_ps = Tf_p(2:nnx+1,1:nny);      
    end
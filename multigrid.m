%function [Tn,res,resc,Xc,Yc,nx,ny] = multigrid(iparam)
clear all
close all
clc

%% Parameters (Taken from mudpack)
% for the finest grid size (nx,ny) on a domain [xa,xb] x [yc,yd]
% ixp -> no.of points of coarsest x grid
% jyq -> no. of points on coarsest y grid
% iex -> exponent used in finding the number of grid points in x
% jey -> exponent used in finding the number of grid points in y
% nfx = ixp*(2**(iex-1)) + 1
% nfy = jyq*(2**(jey-1)) + 1
% ngrid -> number of levels in the multigrid solver
iparam = [0,0,0,0,0,2,2,5,5,0,0];
%
% Boundary condition flags
%
nxa = iparam(2);
nxb = iparam(3);
nyc = iparam(4);
nyd = iparam(5);

%
% Grid size parameters
%
ixp = iparam(6);
jyq = iparam(7);
iex = iparam(8);
jey = iparam(9);
ngrid = max(iex,jey);
nfx = iparam(10);
nfy = iparam(11);

% Set subgrid sizes
for k=1:ngrid
    nxk(k) = ixp*2^(max(k+iex-ngrid,1)-1)+1;
    nyk(k) = jyq*2^(max(k+jey-ngrid,1)-1)+1;
end
nx = nxk(1); %finest grid points in x
ny = nyk(1); %finest grid points in y
psi = ones(nx+2,ny+2);
[x,y,eta,detady,der1,L] = gridCluster(nx,ny);

%% Define surface
phi = surfTop(x,L,nx);

%% Grid clustering function
[Xc,Yc] = meshgrid(x,eta);
for i=1:nx+2
    for j=1:ny
        Yc(j,i) = Yc(j,i) + phi(i);
      %  Y(j,i) = Y(j,i) + phi(i);
    end
end

%% Stencil and Coefficients for the material matrix
[A,mu,psi] = materialMatrix(x,L,nx,ny);
[nw,n,ne,w,ce,e,se,s,sw] = coeff_stencil(A,mu,phi,der1,detady,x,y,nx,ny);

%% Multigrid options - Main steps in Multigrid cycling
% transferc2f - > transfer from coarse mesh to fine
% rst2d -> restrict residual using full weighting operator fine to coarse
% prolong -> interpolate from coarse to fine
% The main idea in multigrid is to perform some sweeps at finest grid,
% restrict the residual to coarse grid, perform some more sweeps at coarser
% level, the interpolate to finer grid and perform some more sweeps (V cycle)

%% Smoothing on the finest grid - Line relaxation in y
T = zeros(nx+2,ny+2); Tn = T; res= T;
niter = 1000;
for iter=1:niter
    for i=2:nx+1
        for j=2:ny
            term_rhs = w(i,j)*T(i-1,j) + e(i,j)*T(i+1,j) + ne(i,j)*T(i+1,j+1)...
                + nw(i,j)*T(i-1,j+1) + se(i,j)*T(i+1,j-1) + sw(i,j)*T(i-1,j-1);
            
            d(j) = -term_rhs;
            
            if (j==2)
                d(j) = d(j) - s(i,j)*T(i,j-1) ;
            elseif (j ==ny)
                d(j) = d(j) + 5; % boundary condition on rhs
            end
        end
        
        Tn = poisson_TDMA(Tn,n,ce,s,i,ny+1,d);
    end
    %bdryCondition(Tn,nx,ny,fBdry);
    Tn(:,1) = 0;
    Tn(1,:) = Tn(nx,:);
    Tn(nx+2,:) = Tn(3,:);
    
    for i=2:nx+1
        for j=2:ny
            res(i,j) = -(nw(i,j)*Tn(i-1,j+1) + n(i,j)*Tn(i,j+1) + ne(i,j)*Tn(i+1,j+1)...
                +   w(i,j)*Tn(i-1,j) + ce(i,j)*Tn(i,j) + e(i,j)*Tn(i+1,j)...
                + sw(i,j)*Tn(i-1,j-1) + s(i,j)*Tn(i,j-1) + se(i,j)*Tn(i+1,j-1));
        end
        res(i,1) = 0; % dirichlet boundary condition
    end
    T = Tn;
end

% %% Restriction operator
% ff = res(2:nx+1,:);
% 
% for l = ngrid:-1:2
%     
%     nnx = nxk(l);
%     nny = nyk(l);
%     
%     % Number of points on coarsest grid
%     ncx = nxk(l-1);
%     ncy = nyk(l-1);
%     resc = zeros(ncx,ncy);
%     
%     % restrict the residual
%     resc = restrct(ff,resc,nnx,nny,ncx,ncy); % (ncx,ncy)
%     res1 = resc;
% 
%     % Regrid and find points on the new grid
%     [xc,yc,etac,detadyc,der1c,L] = gridCluster(ncx,ncy);
%     
%     % Redefine material matrices on new grid
%     [Ac,muc,psic] = materialMatrix(xc,L,ncx,ncy); %(ncx+2,ncy+2)
%     
%     % Coefficient stencil on new grid
%     phic = surfTop(xc,L,ncx);
%     [nwc,nc,nec,wc,cec,ec,sec,sc,swc] = coeff_stencil(Ac,muc,phic,der1c,detadyc,xc,yc,ncx,ncy);
%     
%     % Smooth the residual / Solve for delta-H
%     niter = 10;
%     if (l == 2) niter = 10000; end
%     [deltac] = smoothRes(resc,ncx,ncy,nwc,nc,nec,ec,cec,wc,swc,sc,sec,niter);
%     ff = zeros(ncx,ncy);
%     
%     % Calculate residual for deltaH - recursion here
%     for i=2:ncx+1
%         for j=2:ncy
%             ic = i-1;
%             ip = i+1;
%             im = i-1;
%             if (i==2) im = ncx; end
%             if (i==ncx+1) ip = 3; end
%             ff(ic,j) = resc(ic,j) -(nc(i,j)*deltac(i,j+1) + sc(i,j)*deltac(i,j-1) + cec(i,j)*deltac(i,j)...
%                 + wc(i,j)*deltac(im,j) + ec(i,j)*deltac(ip,j) + nec(i,j)*deltac(ip,j+1)...
%                 + nwc(i,j)*deltac (im,j+1) + sec(i,j)*deltac(ip,j-1) + swc(i,j)*deltac(im,j-1));
%         end
%     end
% end
% delta_cg = deltac(2:ncx+1,1:ncy);
% 
% %% Coarse grid correction of the residual
% % Prolongation operator
% % Consider two meshes VH(fine) and Vh(coarse). The prolongation operator
% % takes vector from coarse to fine grid using linear interpolation i.e.
% %   I(h-H) : VH ---> Vh
% niter = 10;
% for m = ngrid-1%:2
%     nnx = nxk(m+1);
%     nny = nyk(m+1);
%     % Number of points on coarsest grid
%     ncx = nxk(m);
%     ncy = nyk(m);
%     deltaf = zeros(nnx,nny);
%     deltaf = prolong(deltaf,delta_cg,nnx,nny,ncx,ncy);
% 
%     T(2:nnx+1,1:nny) = T(2:nnx+1,1:ny) + deltaf;
%     T(1,1:nny) = T(nnx,1:nny);
%     T(nnx+2,1:ny) = T(3,1:nny);
%     [Tf,rf] = postSmooth(T,nnx,nny,nw,n,ne,e,ce,w,sw,s,se,niter);
% end
% %%
figure
pcolor(Xc(:,2:nx+1),Yc(:,2:nx+1),T(2:nx+1,1:ny)')
%shading interp
% %  figure
% %  surf(Tf(2:nx+1,1:ny))
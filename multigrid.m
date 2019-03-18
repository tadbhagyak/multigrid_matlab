function [Tn,res,Xc,Yc,nx,ny] = multigrid(iparam)
%% Parameters (Taken from mudpack)
% for the finest grid size (nx,ny) on a domain [xa,xb] x [yc,yd]
% ixp -> no.of points of coarsest x grid
% jyq -> no. of points on coarsest y grid
% iex -> exponent used in finding the number of grid points in x
% jey -> exponent used in finding the number of grid points in y
% nfx = ixp*(2**(iex-1)) + 1
% nfy = jyq*(2**(jey-1)) + 1
% ngrid -> number of levels in the multigrid solver

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

L = 0.012;
nx = nxk(ngrid); %finest grid points in x
ny = nyk(ngrid); %finest grid points in y
y = linspace(-2*L,0,ny);
x = linspace(-L,L,nx+2);
[X,Y] = meshgrid(x,y);
dy = abs(y(2) - y(1));
dx = abs(x(2) - x(1));
psi = ones(nx+2,ny+2);

%% Define surface
phi = surfTop(x,L,nx);

%% Grid clustering function
c1y = 0;
yend = 2*L;
eta_2 = (y./yend).^2;
den = 2.0d0 - eta_2;
f = y./(den.^c1y);
der1 = 1./(den.^c1y) + 2*eta_2*c1y./(den.^(c1y+1));
detady = 1./der1;
[Xc,Yc] = meshgrid(x,f);

for i=1:nx+2
    for j=1:ny
        Yc(j,i) = Yc(j,i) + phi(i);
        Y(j,i) = Y(j,i) + phi(i);
    end
end

%% Stencil and Coefficients for the material matrix
A1 = 29176.66;
mu1 = 10670;
mu2 = 2.5;
A2 = 12503.33;

A = A2*ones(nx+2,ny+2);
mu = mu2*ones(nx+2,ny+2);
for i=1:nx+2
    for j=1:ny
        %if ((x(i)^2 + f(j)^2 - 0.5*L^2) <= 0)
        if (and(x(i) > -L/6, x(i) < L/6)) 
            mu(i,j) = mu2;
            A(i,j) = A2;
            psi(i,j) = -1;
        end
    end
end

[nw,n,ne,w,ce,e,se,s,sw] = coeff_stencil(A,mu,phi,der1,detady,x,y,nx,ny);
%% Multigrid options - Main steps in Multigrid cycling
% transferc2f - > transfer from coarse mesh to fine
% rst2d -> restrict residual using full weighting operator fine to coarse
% prolong -> interpolate from coarse to fine
% The main idea in multigrid is to perform some sweeps at finest grid,
% restrict the residual to coarse grid, perform some more sweeps at coarser
% level, the interpolate to finer grid and perform some more sweeps (V cycle)
%% Smoothing on the finest grid - Line relaxation in y
T = zeros(nx+2,ny+2); Tn = T;
niter = 1;
for iter=1:niter
    for i=2:nx+1
        for j=2:ny
            term_rhs = w(i,j)*T(i-1,j) + e(i,j)*T(i+1,j) + ne(i,j)*T(i+1,j+1)...
                + nw(i,j)*T(i-1,j+1) + se(i,j)*T(i+1,j-1) + sw(i,j)*T(i-1,j-1);
            
            d(j) = -term_rhs;
            
            if (j==2)
                d(j) = d(j) - s(i,j)*T(i,j-1) ;
            elseif (j ==ny)
                d(j) = d(j) + 1; % boundary condition on rhs
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
                     + sw(i,j)*Tn(i-1,j+1) + s(i,j)*Tn(i,j-1) + se(i,j)*Tn(i-1,j-1));
        end
    end
    
    T = Tn;
end

%% Restriction operator  
ff =res;
for l = ngrid:-1:ngrid-1
    nnx = nxk(l);
    nny = nyk(l);
    % Number of points on coarsest grid
    ncx = nxk(l-1);
    ncy = nyk(l-1);
    fc = zeros(ncx,ncy);
    fc = restrct(ff,fc,nnx,nny,ncx,ncy);
    ff = fc;
end

%% Coarse grid correction of the residual
%% Prolongation operator
% Consider two meshes VH(fine) and Vh(coarse). The prolongation operator
% takes vector from coarse to fine grid using linear interpolation i.e.
%   I(h-H) : VH ---> Vh
for m = l:ngrid-1
    nnx = nxk(m+1);
    nny = nyk(m+1);
    % Number of points on coarsest grid
    ncx = nxk(m);
    ncy = nyk(m);
    ff = zeros(nnx,nny);
    ff = prolong(ff,fc,nnx,nny,ncx,ncy);
    fc = ff;
end
    

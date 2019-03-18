clear all
close all
clc

%% ------------------------------------------------------------------------
%   Multigrid for poisson problem
%   - V cycles, W cycles and Full multigrid cycles
%   - Line relxation in y for sweeping
%   - 4 level of grids
%   - Poisson solver for velocity 
%------------------------------------------------------------------------
% for the finest grid size (nx,ny) on a domain [xa,xb] x [yc,yd]
% ixp -> no.of points of coarsest x grid
% jyq -> no. of points on coarsest y grid
% iex -> exponent used in finding the number of grid points in x
% jey -> exponent used in finding the number of grid points in y
% nfx = ixp*(2**(iex-1)) + 1
% nfy = jyq*(2**(jey-1)) + 1
% ngrid -> number of levels in the multigrid solver

iparam = [0,0,0,0,0,2,2,5,5,0,0];
[f,res,resc,X,Y,nx,ny] = multigrid(iparam);

%%Plotting and post processing
%%
figure
pcolor(X(:,2:nx+1),Y(:,2:nx+1),f(2:nx+1,1:ny)')
shading interp
colormap jet
colorbar
grid on
%hold on
%contour(Xc(:,2:nx+1),Yc(:,2:nx+1),psi(2:nx+1,:)',[-0.001 0.001],'w')

%%
% figure
% fbar = f + phi(nx/2);
% Taly = 10*(fbar - fbar(1))./(fbar(end) - fbar(1));
% yctr = Yc(:,nx/2);
% plot(fbar,Taly,'LineWidth',2)
% hold on
% plot(yctr,T(nx/2,:),'-o')

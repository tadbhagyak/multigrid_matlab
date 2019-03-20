function [x,y,eta,detady,der1,L] = gridCluster(nx,ny)
L = 0.012;
y = linspace(-2*L,0,ny);
x = linspace(-L,L,nx+2);
[X,Y] = meshgrid(x,y);
dy = abs(y(2) - y(1));
dx = abs(x(2) - x(1));
%% Grid clustering function
c1y = 0;
yend = 2*L;
eta_2 = (y./yend).^2;
den = 2.0d0 - eta_2;
eta = y./(den.^c1y);
der1 = 1./(den.^c1y) + 2*eta_2*c1y./(den.^(c1y+1));
detady = 1./der1;
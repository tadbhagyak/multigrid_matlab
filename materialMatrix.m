function [A,mu,psi] = materialMatrix(x,L,nx,ny)

A1 = 100;%29176.66;
mu1 = 100;%10670;
mu2 = 100;%2.5;
A2 = 100;%12503.33;

A = A1*ones(nx+2,ny);
mu = mu1*ones(nx+2,ny);
psi = ones(nx+2,ny);
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
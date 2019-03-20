function [A,mu,psi] = materialMatrix(x,L,nx,ny)

A1 = 29176.66;
mu1 = 10670;
mu2 = 2.5;
A2 = 12503.33;

A = A1*ones(nx+2,ny);
mu = mu1*ones(nx+2,ny);
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
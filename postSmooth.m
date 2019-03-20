function [Tn,res] = postSmooth(T,nx,ny,nw,n,ne,e,ce,w,sw,s,se,niter)
Tn = T;
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
                + se(i,j)*Tn(i+1,j-1) + s(i,j)*Tn(i,j-1) + sw(i,j)*Tn(i-1,j-1));
        end
        res(i,1) = 0; % dirichlet boundary condition
    end
    T = Tn;
end
function [deltac] = smoothRes(resc,ncx,ncy,nwc,nc,nec,ec,cec,wc,swc,sc,sec,niter)

deltac = zeros(ncx+2,ncy+2);
rescn = deltac;
for iter=1:niter
    for i=2:ncx+1
        for j=2:ncy
            ic = i-1;
            ip = i+1;
            im = i-1;
            if (i==2) im = ncx; end
            if (i==ncx+1) ip= 3; end
            term_rhs = wc(i,j)*deltac(im,j) + ec(i,j)*deltac(ip,j) + nec(i,j)*deltac(ip,j+1)...
                + nwc(i,j)*deltac (im,j+1) + sec(i,j)*deltac(ip,j-1) + swc(i,j)*deltac(im,j-1);
            dc(j) = resc(ic,j) - term_rhs;
            if (j==2)
                dc(j) = dc(j) - sc(i,j)*deltac(i,j-1) ;
            elseif (j ==ncy)
                %    dc(j) = dc(j) - nc(i,j)*deltac(i,j+1) ;
            end
        end
        rescn = poisson_TDMA(rescn,nc,cec,sc,i,ncy+1,dc);
    end
    deltac = rescn;
end

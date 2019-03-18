function [nwc,nc,nec,wc,cec,ec,sec,sc,swc] = coeff_stencil_coarse(nw,n,ne,w,ce,e,se,s,sw,nxc,nyc,nx,ny)
nwc = zeros(nxc,nyc);nc=nwc;nec=nwc;wc=nwc;cec=nwc;ec=nwc;sec=nwc;
sc=nwc; swc=nwc;

for ic=1:nxc
    for jc=1:nyc
        i = 2*ic-1;
        j=  2*jc-1;
        if (i > nx+2) 
            break;
        end
        if (j > ny) 
            break; 
        end;
            
        nwc(ic,jc) = nw(i,j);
        nc(ic,jc) = n(i,j);
        nec(ic,jc) = ne(i,j);
        wc(ic,jc) = w(i,j);
        cec(ic,jc) = ce(i,j);
        ec(ic,jc) = e(i,j);
        sec(ic,jc) = se(i,j);
        sc(ic,jc) = s(i,j);
        swc(ic,jc) = sw(i,j);
    end
end
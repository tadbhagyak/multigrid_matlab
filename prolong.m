function [ff] = prolong(ff,fc,nnx,nny,ncx,ncy)
for i=1:ncx
    for j=1:ncy
        ff(2*i-1,2*j-1) = fc(i,j);
        
        if (i<ncx)
            ff(2*i,2*j-1) = 0.5*(fc(i,j) + fc(i+1,j));
        end
        if (j<ncy)    
            ff(2*i-1,2*j) = 0.5*(fc(i,j) + fc(i,j+1));
        end
        if (and(i< ncx, j <ncy))
            if (2*i > nnx) break; end
            if (2*j > nny) break; end
            ff(2*i,2*j) = 0.25*(fc(i,j) + fc(i+1,j) + fc(i,j+1) + fc(i+1,j+1));
        end
    end
end

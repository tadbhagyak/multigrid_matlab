function [nw,n,ne,w,ce,e,se,s,sw] = coeff_stencil(A,mu,phi,der1,detady,x,y,nx,ny)

nw = zeros(nx+2,ny);
n=nw;ne=nw;w=nw;ce=nw;e=nw;se=nw;s=nw;sw=w;

% Formulate the stencil matrix
for i=2:nx+1
    for j=2:ny
        dx = abs(x(2)-x(1));
        fx1 = (phi(i) - phi(i-1))/dx;
        fx2 = (phi(i+1) - phi(i))/dx;
        fx3 = fx1;
        fx4 = fx2;
        
        dy1 = 0;
        dy2 = 0;
        if (j > 1)  dy1 = abs(y(j) - y(j-1));end
        if (j < ny) dy2 = abs(y(j+1) - y(j));end
        
        % Coefficents associated with phi_x T_X
        c1 = der1(j)*dy1/(4*dx) + der1(j-1)*dy1/(12*dx);
        ce_1 = der1(j)*dy1/(12*dx) + der1(j-1)*dy1/(12*dx);
        if (j<ny)
            c5 = der1(j)*dy2/(4*dx) + der1(j+1)*dy2/(12*dx);
            ce_5 = der1(j)*dy2/(12*dx) + der1(j+1)*dy2/(12*dx);
        end
        
        % Coefficents associated with phi_eta T_eta
        detady1 = (detady(j-1) + detady(j))/(2);
        if (j<ny) detady2 = (detady(j) + detady(j+1))/(2);end
        
        m1 = mu(i-1,j-1);
        m2 = mu(i,j-1);
        m3 = mu(i-1,j);
        m4 = mu(i,j);
        
        c2_a1 = dx/(6*dy1)*detady1*(A(i-1,j-1) + mu(i-1,j-1)*fx1^2);
        c2_a2 = dx/(6*dy1)*detady1*(A(i,j-1) + mu(i,j-1)*fx2^2);
        c2_a3 = dx/(6*dy2)*detady2*(A(i-1,j) + mu(i-1,j)*fx3^2);
        c2_a4 = dx/(6*dy2)*detady2*(A(i,j) + mu(i,j)*fx4^2);
        
        sw(i,j) = -m1*ce_1 - c2_a1 + 2*m1*fx1/4;
        se(i,j) = -m2*ce_1 - c2_a2 - 2*m2*fx2/4;
        nw(i,j) = -m3*ce_5 - c2_a3 - 2*m3*fx3/4;
        ne(i,j) = -m4*ce_5 - c2_a4 + 2*m4*fx4/4;
        
        w(i,j) = -m1*c1 - m3*c5 + c2_a1 + c2_a3;
        e(i,j) = -m2*c1 - m4*c5 + c2_a2 + c2_a4;
        
        s(i,j) = m1*ce_1 + m2*ce_1 - 2*c2_a1 - 2*c2_a2;
        n(i,j) = m3*ce_5 + m4*ce_5 - 2*c2_a3 - 2*c2_a4;
        
        ce(i,j) = (m1*c1 + m2*c1 + m3*c5 + m4*c5 + 2*c2_a1 + 2*c2_a2...
            + 2*c2_a3 + 2*c2_a4) - 2*(-fx1*m1 + fx2*m2 + fx3*m3 - fx4*m4)/4;
        
        if (j==ny)
            sw(i,j) = -m1*ce_1 - c2_a1 + 2*m1*fx1/4;
            se(i,j) = -m2*ce_1 - c2_a2 - 2*m2*fx2/4;
            nw(i,j) = 0;
            ne(i,j) = 0;
            w(i,j) = -m1*c1 + c2_a1;
            e(i,j) = -m2*c1 + c2_a2;
            s(i,j) = m1*ce_1 + m2*ce_1 - 2*c2_a1 - 2*c2_a2;
            n(i,j) = 0;
            ce(i,j) = (m1*c1 + m2*c1 + 2*c2_a1 + 2*c2_a2) - 2*(-fx1*m1 + fx2*m2)/4;
        end
        
        if (j==1)
            sw(i,j) = 0;%-m1*ce_1 - c2_a1 + 2*m1*fx1/4;
            se(i,j) = 0;%-m2*ce_1 - c2_a2 - 2*m2*fx2/4;
            nw(i,j) = -m3*ce_5 - c2_a3 - 2*m3*fx3/4;
            ne(i,j) = -m4*ce_5 - c2_a4 + 2*m4*fx4/4;
            
            w(i,j) =  - m3*c5  + c2_a3;
            e(i,j) =  - m4*c5 + c2_a4;
            
            s(i,j) = 0;%m1*ce_1 + m3*ce_1 - 2*c2_a1 - 2*c2_a2;
            n(i,j) = m3*ce_5 + m4*ce_5 - 2*c2_a3 - 2*c2_a4;
            
            ce(i,j) = (m3*c5 + m4*c5 ...
                + 2*c2_a3 + 2*c2_a4) - 2*(fx3*m3 - fx4*m4)/4;
        end
    end
end
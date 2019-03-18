function [fn] = poisson_TDMA(fn,c,b,a,i,ny,d)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%modified coefficients
c1(2) = c(i,2)/b(i,2);
d1(2) = d(2)/b(i,2);

for j = 3:ny-1
    c1(j) = c(i,j)/( b(i,j) - a(i,j)*c1(j-1));
    d1(j) = (d(j) - a(i,j)*d1(j-1))/( b(i,j) - a(i,j)*c1(j-1));
    
end

fn(i,ny-1) = d1(ny-1);

%==============================Back Substitution==============================

for j = ny-2:-1:2
    fn(i,j) = -c1(j)*fn(i,j+1) + d1(j);
end

end


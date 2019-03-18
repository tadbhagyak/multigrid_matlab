function [phi] = surfTop(x,L,nx)
%% Surface function
phi = zeros(nx+2,1);
amp = -0.002;
for i =2:nx+1
    % phi(i) = amp*exp(-66.0*((x(i)-0.002)/L)^2)...
    %  +amp*exp(-66.0*((x(i)+0.002)/L)^2);
    phi(i) = amp*exp(-10.0*((x(i))/L)^2);
end
phi(1) = phi(nx);
phi(nx+2) = phi(3);
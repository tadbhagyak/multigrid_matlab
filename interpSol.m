function [Tf] = interpSol(Tc,nxk,nyk,kcycle)

nnx = nxk(kcycle+1);
nny = nyk(kcycle+1);
ncx = nxk(kcycle);
ncy = nyk(kcycle);
Tf = zeros(nnx,nny);
Tf = prolong(Tf,Tc,nnx,nny,ncx,ncy);
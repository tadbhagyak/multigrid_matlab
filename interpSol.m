function [Tf] = interpSol(Tc,nxk,nyk,cur_grid)

nnx = nxk(cur_grid+1);
nny = nyk(cur_grid+1);
ncx = nxk(cur_grid);
ncy = nyk(cur_grid);
Tf = zeros(nnx,nny);
Tf = prolong(Tf,Tc,nnx,nny,ncx,ncy);
% Function to compute the vd and vq in one extreme according to p and q sent in the other extreme. 
% It receives the voltage magnitude of the other extreme (considering this
% as a reference Vd = vmag, Vq =0

function [vd, vq, vmag] = pq2vdq(r, x, v, p, q)

vd = v*(1 + r*p/v + x*q/v);
vq = x*p/v - r*q/v;
vmag = (vd^2 + vq^2)^0.5;
return
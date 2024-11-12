% Function that calculates the active and reactive power of a 
% three-phase RL load (R-L in series) given an active and reactive power at nominal voltage.
% p: active power [W]
% q: reactive power [VAR], if positive then inductive, if negative then
% capacitive
% v: nominal voltage [V]

function [p,q] = z2pq(r,x, v)

p = v^2 * r/(r^2 + x^2);
q = v^2 * x/(r^2 + x^2);

end
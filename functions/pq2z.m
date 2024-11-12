% Function that calculates the resistance (r) and reactance (x) of a
% three-phase RL load (R-L in series) given an active and reactive power at nominal voltage.
% p: active power [W]
% q: reactive power [VAR], if positive then inductive, if negative then
% capacitive
% v: nominal voltage [V]

function [r,x] = pq2z(p,q, v)
if p>0
    ang = atan(q/p);
    z = v^2*cos(ang)/p;
    r = z*cos(ang);
    x = z*sin(ang);
else
    r = 0;
    x = v^2/q ;
end
end
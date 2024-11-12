function [ssr, ssrt, Rl_ohm, Rc_ohm, L_H, C_F] = RLseries_RCsh(Load, Grid)

[Rl_ohm, XL]  = pq2z(Load.PL, Load.QL, Grid.Vbase);
[~, XC]  = pq2z(0, Load.QC, Grid.Vbase);
[Rc_ohm, ~]  = pq2z(Load.PC, 0, Grid.Vbase);

L_H        = XL/Grid.w ;
C_F        = -1/(Grid.w*XC);

Rl     = Rl_ohm/Grid.Zbase;
Rc     = Rc_ohm/Grid.Zbase;
L      = L_H*Grid.w/Grid.Zbase;
C      = C_F*Grid.w*Grid.Zbase;

x = ["A"; "B";  "C"; "D"];
y = ["Ath"; "Bth";  "Cth"; "Dth"];


% State-space realization in p.u. (as it is)
syms rl lc rc c w wo

ssr.Ath = wo*[  -rl/lc    w      1/lc           0; 
                -w      -rl/lc   0              1/lc;
                -1/c    0        -1/(rc*c)        w;
                0      -1/c      -w           -1/(rc*c)];

ssr.Bth = wo*[ 0 0; 0 0; 1/c 0;0 1/c];
ssr.Cth = eye(4,4);
ssr.Dth = zeros(4,2);

for j=1:length(x)
    ssr.(x(j)) = double(subs(ssr.(y(j)), ...
              [rl,       lc,    rc,   c, w, wo], ...
              [Rl,       L,     Rc,   C, 1, Grid.w]));
end

% State-space realization in time-dimensionless
syms rl lc rc c w

ssrt.Ath =    [  -rl/lc    w      1/lc           0; 
                -w      -rl/lc   0              1/lc;
                -1/c    0        -1/(rc*c)        w;
                0      -1/c      -w           -1/(rc*c)];

ssrt.Bth = [ 0 0; 0 0; 1/c 0;0 1/c];
ssrt.Cth = eye(4,4);
ssrt.Dth = zeros(4,2);

for j=1:length(x)
    ssrt.(x(j)) = double(subs(ssrt.(y(j)), ...
              [rl,       lc,    rc,   c, w, wo], ...
              [Rl,       L,     Rc,   C, 1, Grid.w]));
end


end

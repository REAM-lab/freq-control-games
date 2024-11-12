function [ssr, ssrt, R_ohm, L_H, C_F, R, L, C] = RLseries_Csh(Load, Grid)

[R_ohm, XL]  = pq2z(Load.P, Load.QL, Grid.Vbase);
[~, XC]  = pq2z(0, Load.QC, Grid.Vbase);

L_H        = XL/Grid.w ;
C_F        = -1/(Grid.w*XC);

R     = R_ohm/Grid.Zbase;
L     = L_H*Grid.w/Grid.Zbase;
C     = C_F*Grid.w*Grid.Zbase;

x = ["A"; "B";  "C"; "D"];
y = ["Ath"; "Bth";  "Cth"; "Dth"];


% State-space realization in p.u. (as it is)
syms r lc c w wo

ssr.Ath = wo*[  -r/lc    w      1/lc      0; 
                -w      -r/lc   0        1/lc;
                -1/c    0      0        w;
                0      -1/c    -w       0];

ssr.Bth = wo*[ 0 0; 0 0; 1/c 0;0 1/c];
ssr.Cth = eye(4,4);
ssr.Dth = zeros(4,2);

for j=1:length(x)
    ssr.(x(j)) = double(subs(ssr.(y(j)), ...
              [r,       lc,      c, w, wo], ...
              [R,       L,      C, 1, Grid.w]));
end

% State-space realization in time-dimensionless
syms r lc c w

ssrt.Ath = [  -r/lc    w      1/lc      0; 
                -w      -r/lc   0        1/lc;
                -1/c    0      0        w;
                0      -1/c    -w       0];

ssrt.Bth = [ 0 0; 0 0; 1/c 0;0 1/c];
ssrt.Cth = eye(4,4);
ssrt.Dth = zeros(4,2);

for j=1:length(x)
    ssrt.(x(j)) = double(subs(ssrt.(y(j)), ...
              [r, lc,  c, w], ...
              [R, L,  C, 1]));
end


end

function [ssr, ssrt, R_ohm, L_H]  = RLbranch(TX, Grid)

x = ["A"; "B";  "C"; "D"];
y = ["Ath"; "Bth";  "Cth"; "Dth"];

syms R L w wo

ssr.Ath = wo*[-R/L w;-w -R/L];
ssr.Bth = wo*[1/L 0 -1/L 0;0 1/L 0 -1/L];
ssr.Cth = eye(2,2);
ssr.Dth = zeros(2,4);

for j=1:length(x)
    ssr.(x(j)) = double(subs(ssr.(y(j)), ...
              [R,       L,  w, wo], ...
              [TX.R,    TX.L, 1, Grid.w]));
end

syms R L w

ssrt.Ath = [-R/L w;-w -R/L];
ssrt.Bth = [1/L 0 -1/L 0;0 1/L 0 -1/L];
ssrt.Cth = eye(2,2);
ssrt.Dth = zeros(2,4);

for j=1:length(x)
    ssrt.(x(j)) = double(subs(ssrt.(y(j)), ...
              [R,       L,  w], ...
              [TX.R,    TX.L, 1]));
end

R_ohm = TX.R * Grid.Zbase ; 
L_H = TX.L * Grid.Zbase/Grid.w;

end
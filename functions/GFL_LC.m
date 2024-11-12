% University of California San Diego
% Jacobs School of Engineering
% Authors: Paul Serna-Torre, ..... (looking for someone who can contribute :-) )

% Function that provides the state-space realization of the voltage control
% scheme of a grid-following IBR with dq-frame controllers
% IBR is a structure that has the parameters in p.u.

% Three state-space realization are provided by the function:
%   1) ssr : it is the state-space realization with units in p.u.
%   2) ssrt: it is the state-space realization with units in p.u. and
%            time-dimensionless, where t_new = t * wo. This change of
%            variable is effective to reduce the large number variations in
%            matrix A. The eigenvalues change in magnitude with respect to
%            1) but the damping does not change and the frequency natural
%            is rescale by wo.
%   3) ssrT: it is the state-space realization that takes the ssr of 2) and
%            apply a transformation x_new = T*x, where T is the transformation
%            matrix. This transformation is necessary to better the condition of the
%            matrix A. In 1) and 2) the matrix A may be ill-conditioned which can
%            mess the calculation of rank and controllability. The matrix T
%            is given by the user, and must saved under the name TrM as a
%            field in the structure IBR.
%
%   In 1), 2) and 3), theoretical matrices are also given. The user can access 
%   by calling, for example, IBR.ssr.Ath.
%


function [ssr, ssrt, ssrT] = GFL_LC(IBR)

% u = [iref_dq, v_dq], y = i_dq, x = [phi_dq, e_dq, i_dq]


% Nomenclature
% w     : frequency [p.u.] 
% wo    : nominal frequency [Hz]
% Kp    : Gain of proportional action in PI of voltage loop
% Ki    : Gain of integral action in PI of voltage loop
% Rc    : Resistance 1 of LCL filter [p.u.]
% Lc    : Inductance 1 of LCL filter [p.u.]
% Tff   : Time constant of first-order filter for feed-forward compensation

% Note: Though the current function does not care of the Ohm, H, or F, we
%       still provide a guidance of how they should be calculated outside of this
%       function. We prefer to use p.u. since the values of R, L, and C are
%       normalized.
%       Rf[p.u.] = Rf[Ohm]/Zbase,
%       Lf[p.u.] = Lf[H]*wo/Zbase,
%       Cf[p.u.] = Cf[F]/C[base] = Cf[F]*wo*Zbase ,
%           Notice that while Rf and Lf are in p.u. of Zbase, Cf is ...
%           in p.u. of 1/Zbase. So Cf[p.u] is usually a small value like
%           0.001;
%           

% The state-space realization was derived using the CCM method and joining
% the ssr of the LC filter and controllers

x = ["A"; "B";  "C"; "D"];
y = ["Ath"; "Bth";  "Cth"; "Dth"];

%% State-space realization with p.u. (as it is)
syms w wo Kp Ki Rc Lc Tff

% Current controller
CC.ssr.Ath = zeros(2,2);
CC.ssr.Bth = [1 0;0 1];
CC.ssr.Cth = [Ki 0;0 Ki];
CC.ssr.Dth = [Kp 0;0 Kp];

% Feed-forward compensation
FF.ssr.Ath = -1/Tff * eye(2,2);
FF.ssr.Bth = 1/Tff*eye(2,2);
FF.ssr.Cth = eye(2,2);
FF.ssr.Dth = zeros(2,2);

% LC filter 
LC.ssr.Ath = wo*[-Rc/Lc  w ;-w  -Rc/Lc];
LC.ssr.Bth = wo*[1/Lc 0 -1/Lc 0;0 1/Lc 0 -1/Lc];
LC.ssr.Cth = eye(2,2);
LC.ssr.Dth = zeros(2,4);


R_3 =  [  zeros(1,4) -1 0;
            zeros(1,5) -1;
            zeros(2,6);
            1 0 1 0 0 -w*Lc;
            0 1 0 1 w*Lc 0;
            zeros(2,6)];
R_2 = [eye(4,4);
            zeros(2,4);
            zeros(2,2) eye(2,2)];
R_1 = [zeros(2,4) eye(2,2)];
R_0 = zeros(2,4);

% Theoretical [Ath,Bth,Cth,Dth] of the IBR
ssr = ccm_ssr(R_3, R_2, R_1, R_0, y, y, CC.ssr, FF.ssr, LC.ssr);

% State-space realization with numerical values
% In the following, numerical values are provided to the symbolic variables. Everything is
% p.u., except wo which is 60 Hz or 50 Hz. It is assumed that the
% frequency is fixed at nominal value. Otherwise, the state-space
% realization above would need to consider the frequency w as a state.

for j=1:length(x)
    ssr.(x(j)) = double(subs(ssr.(y(j)), ...
              [w, wo,    Rc,     Lc,     Kp,     Ki,     Tff], ...
              [1, IBR.w, IBR.Rc, IBR.Lc, IBR.Kp, IBR.Ki, IBR.Tf]));
end

%% State-space realization with p.u. in time-dimensionless
syms w wo Kp Ki Rc Lc Tff

% Current controller
CCt.ssr.Ath = zeros(2,2);
CCt.ssr.Bth = [1/wo 0;0 1/wo];
CCt.ssr.Cth = [Ki 0;0 Ki];
CCt.ssr.Dth = [Kp 0;0 Kp];

% Feed-forward compensation
FFt.ssr.Ath = -1/(wo*Tff) * eye(2,2);
FFt.ssr.Bth = 1/(wo*Tff)*eye(2,2);
FFt.ssr.Cth = eye(2,2);
FFt.ssr.Dth = zeros(2,2);

% LC filter 
LCt.ssr.Ath = [-Rc/Lc  w ;-w  -Rc/Lc];
LCt.ssr.Bth = [1/Lc 0 -1/Lc 0;0 1/Lc 0 -1/Lc];
LCt.ssr.Cth = eye(2,2);
LCt.ssr.Dth = zeros(2,4);

ssrt = ccm_ssr(R_3, R_2, R_1, R_0, y, y, CCt.ssr, FFt.ssr, LCt.ssr);

% Numerical values are provided to the symbolic variables. Everything is
% p.u., except the wo which is 60 Hz or 50 Hz. It is assumed that the
% frequency is fixed at nominal value. Otherwise, the state-space
% realization above would need to consider the frequency w as a state.

for j=1:length(x)
    ssrt.(x(j)) = double(subs(ssrt.(y(j)), ...
              [w, wo,    Rc,     Lc,     Kp,     Ki,     Tff], ...
              [1, IBR.w, IBR.Rc, IBR.Lc, IBR.Kp, IBR.Ki, IBR.Tf]));
end

%% State-space realization with transformation 

if isfield(IBR,'TrM')
    
    ssrT.Ath = IBR.TrM * ssrt.Ath * (IBR.TrM)^(-1);
    ssrT.Bth = IBR.TrM * ssrt.Bth;
    ssrT.Cth = ssrt.Cth * (IBR.TrM)^(-1);
    ssrT.Dth = ssrt.Dth;
    
    for j=1:length(x)
        ssrT.(x(j)) = double(subs(ssrT.(y(j)), ...
              [w, wo,    Rc,     Lc,     Kp,     Ki,     Tff], ...
              [1, IBR.w, IBR.Rc, IBR.Lc, IBR.Kp, IBR.Ki, IBR.Tf]));
    end


else
    ssrT='No Transformation Matrix (TrM) exists in structure IBR';

end

end



clear 
clc

%% Info: 
% The file contains the state-space realization of 1 grid-following IBR, 1
% grid-forming IBR, 1 RLC load, 1 TX branch connecting the system with the
% upper-level transmission system.

% Notes:
% (1) If transformer is not used, then set R1, X1, ... to zero. Do
% and remove it from the simulink model. Do not set the Snom of the TXR to
% zero since it is a denominator for calculation Rc.

% (2) The Internal impendance (Rci and Lci) are joined with the impedance
% of the transformer to become the coupling impedance (Rc and Lc). It
% facilitates the calculation of the SSR.

clear 
clc

%% Grid parameters
Grid.Sbase = 100e6;
Grid.Vbase = 230e3;
Grid.Ibase = Grid.Sbase/(sqrt(3)*Grid.Vbase);
Grid.Zbase = Grid.Vbase/(sqrt(3)*Grid.Ibase);
Grid.f = 60;
Grid.w = 2*pi*Grid.f;
Ts = 5e-6;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SG1 (BUS 4)

SG(1).Sbase = Grid.Sbase ;
SG(1).Vbase = 13.8e3; 
SG(1).H     = 4;
SG(1).KD    = 1;

SG(1).TXR.Snom  = Grid.Sbase;
SG(1).TXR.R1    = 0.001;
SG(1).TXR.L1    = 0.0576/2;
SG(1).TXR.R2    = 0.001;
SG(1).TXR.L2    = 0.0576/2;
SG(1).Rp        = 0.05;


%% IBR-GFL1 (BUS 4)

IBR(1) = struct(...
        'Sbase'  , 100e6,...                                % Three-phase base power [VA]
        'Vbase'  , 480,...                                % Nominal line-line rms voltage [V]
        'w'      , Grid.w, ...                            % Nominal frequency [rad/s]
        'P_set'  , 0 ,...                                 % Active power setpoint [p.u.]
        'Rp'     , 0 ,...                                 % Droop 
        'Vdc'    , 4000,...                               % DC voltage [V]
        'Rf'     , 0.01,...                               % Resistance 1 of LCL filter [p.u.]
        'Lf'     , 0.25,...                                % Inductance 1 of LCL filter [p.u.]
        'Cf'     , 0.9,...                                % Capacitance of LCL filter [p.u.]
        'Rci'    , 0,...                                  % Internal resistance 2 of LCL filter [p.u.]
        'Lci'    , 0,...                                  % Internal inductance 2 of LCL filter [p.u.]
        'wc'     , 31.41,...                              % Time constant of filter [s]
        'FF'     , 0.75,...                               % Gain of feed-forward compensation
        'Kpv'    , 0.8,...                                % Gain of proportional controller in voltage loop
        'Kiv'    , 150,...                                % Gain of integral controller in voltage loop
        'Kpc'    , 0.8,...                                % Gain of proportional controller in current loop
        'Kic'    , 150,...                                % Gain of integral controller in current loop
        'Tf'    , 10^(-3), ... %Not used
        'Kp'    , 8.0, ... %Not used
        'Ki'    , 25,... %Not used
        'TrM'    , blkdiag(Grid.w*eye(4,4), 1*eye(2,2)),... % Transformation for time-dimensionless state-space realization
        ...
        TXR = struct(...
        'Snom'   ,   100e6,...         % Transformer nominal power [VA]
        'R1'     ,   0.001,...          % Primary-side R1 [p.u. in base of TXR]
        'L1'     ,   0.0576/2,...          % Primary-side X1 [p.u. in base of TXR]
        'R2'     ,   0.001,...          % Primary-side R2 [p.u. in base of TXR]
        'L2'     ,   0.0576/2));           % Primary-side X2 [p.u. in base of TXR]

% Extra calculations
IBR(1).Ibase = IBR(1).Sbase/(sqrt(3)*IBR(1).Vbase);    % Nominal line current [A]
IBR(1).Zbase = IBR(1).Vbase/(sqrt(3)*IBR(1).Ibase);    % Base impendance [Ohm]
IBR(1).Rc     = ((IBR(1).TXR.R1 + IBR(1).TXR.R2)*IBR(1).Sbase/IBR(1).TXR.Snom) ...
                + IBR(1).Rf;                          % Resistance 2 of LCL filter [p.u.]
IBR(1).Lc     = ((IBR(1).TXR.L1 + IBR(1).TXR.L2)*IBR(1).Sbase/IBR(1).TXR.Snom) ...
                + IBR(1).Lf;                          % Resistance 2 of LCL filter [p.u.]

[IBR(1).ssr, IBR(1).ssrt, IBR(1).ssrT] = GFL_LC(IBR(1)); % Compute state-space realizations
IBR(1).ssr.B = IBR(1).ssr.B(:,1:2);
IBR(1).ssr.C = [zeros(1,4) 1 0];  IBR(1).ssrt.C = IBR(1).ssr.C ; 
IBR(1).ssr.D = zeros(2,4); IBR(1).ssrt.D = zeros(1,2); 


%% IBR-GFL2 (BUS 4)

IBR(2) = IBR(1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TX1 (BUS 4 - BUS 9)

TX(1).R = 0.01;
TX(1).L = 0.085;
TX(1).Y = TX(1).L/(TX(1).R^2 + TX(1).L^2);

[TX(1).ssr, TX(1).ssrt, TX(1).R_ohm, TX(1).L_H]  = RLbranch(TX(1), Grid);

%% SG2 (BUS 9)

SG(2) = SG(1);
SG(2).H = 30/100*2;
SG(2).Rp = 0.05*100/30;

%% LOAD (BUS 9)

L(1) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 160e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 60e6);       % Reactive power [VA])

%% LOAD (BUS 9)

L(2) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 20e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 0e6);       % Reactive power [VA])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TX2 (BUS 9 - BUS 8)

TX(2).R = 0.032;
TX(2).L = 0.161;
TX(2).Y = TX(2).L/(TX(2).R^2 + TX(2).L^2);

[TX(2).ssr, TX(2).ssrt, TX(2).R_ohm, TX(2).L_H]  = RLbranch(TX(2), Grid);

%% SG3 (BUS 8)

SG(3) = SG(1);

%% IBR-GFL2 (BUS 8)

IBR(2) = IBR(1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TX3 (BUS 8 - BUS 7)

TX(3).R = 0.0085;
TX(3).L = 0.072;
TX(3).Y = TX(3).L/(TX(3).R^2 + TX(3).L^2);

[TX(3).ssr, TX(3).ssrt, TX(3).R_ohm, TX(3).L_H]  = RLbranch(TX(3), Grid);

%% SG4 (BUS 7)

SG(4) = SG(1);

%% LOAD (BUS 7)

L(3) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 50e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 0e6);       % Reactive power [VA])

%% LOAD (BUS 7)

L(4) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 50e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 0e6);       % Reactive power [VA])

%% IBR-GFL4 (BUS 8)

IBR(4) = IBR(1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TX4 (BUS 7 - BUS 6)

TX(4).R = 0.0119;
TX(4).L = 0.1008;
TX(4).Y = TX(4).L/(TX(4).R^2 + TX(4).L^2);

[TX(4).ssr, TX(4).ssrt, TX(4).R_ohm, TX(4).L_H]  = RLbranch(TX(4), Grid);

%% SG5 (BUS 6)

SG(5) = SG(1);

%% IBR-GFL3 (BUS 6)

IBR(3) = IBR(1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUS 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TX5 (BUS 6 - BUS 5)

TX(5).R = 0.039;
TX(5).L = 0.17;
TX(5).Y = TX(5).L/(TX(5).R^2 + TX(5).L^2);

[TX(5).ssr, TX(5).ssrt, TX(5).R_ohm, TX(5).L_H]  = RLbranch(TX(5), Grid);

%% SG6 (BUS 5)

SG(6) = SG(1);
SG(6).H = 30/100*2;
SG(6).Rp = 0.05*100/30;

%% LOAD (BUS 5)

L(5) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 150e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 20e6);       % Reactive power [VA])

%% LOAD (BUS 5)

L(6) = struct(...
        'Sbase' , Grid.Sbase,...  % Nominal 3ph power [VA]
        'Vbase' , Grid.Vbase,...  % Nominal line-line rms voltage [V]
        'P'     , 50e6,...       % Active power [W]
        'QL'    , 0e6,...       % Reactive power [VA] 
        'QC'    , 0e6);       % Reactive power [VA])

%% TX6 (BUS 5 - BUS 4)

TX(6).R = 0.017;
TX(6).L = 0.092;
TX(6).Y = TX(6).L/(TX(6).R^2 + TX(6).L^2);

[TX(6).ssr, TX(6).ssrt, TX(6).R_ohm, TX(6).L_H]  = RLbranch(TX(6), Grid);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nibrs   = 4 ; %GLFs
nbus    = 6;
ns_ibr  = size(IBR(1).ssr.A,1);
nu_ibr  = size(IBR(1).ssr.B,2);

KD_m = [];
for i=1:nbus
    KD_m = blkdiag(KD_m, -(SG(i).KD+1/SG(i).Rp)/(2*SG(i).H)) ;
end
        %4                                  
L_m = [-1/(2*SG(1).H)*(TX(1).Y + TX(6).Y)     1/(2*SG(1).H)*TX(1).Y                     0                                   0                                   0                                   1/(2*SG(1).H)*TX(6).Y;
        1/(2*SG(2).H)*TX(1).Y                 -1/(2*SG(2).H)*(TX(1).Y + TX(2).Y)        1/(2*SG(2).H)*TX(2).Y               0                                   0                                   0;
        0                                     1/(2*SG(3).H)*TX(2).Y                     -1/(2*SG(3).H)*(TX(2).Y + TX(3).Y)  1/(2*SG(3).H)*TX(3).Y               0                                   0;
        0                                     0                                         1/(2*SG(4).H)*TX(3).Y               -1/(2*SG(4).H)*(TX(3).Y + TX(4).Y)  1/(2*SG(4).H)*TX(4).Y               0;
        0                                     0                                         0                                   1/(2*SG(5).H)*TX(4).Y               -1/(2*SG(5).H)*(TX(4).Y + TX(5).Y)  1/(2*SG(5).H)*TX(5).Y;
        1/(2*SG(6).H)*TX(6).Y                 0                                         0                                   0                                   1/(2*SG(6).H)*TX(5).Y               -1/(2*SG(6).H)*(TX(5).Y + TX(6).Y)];

IBR_m = [ [1/(2*SG(1).H)*IBR(1).ssr.C ; zeros(5,ns_ibr)] [zeros(2,ns_ibr); 1/(2*SG(3).H)*IBR(2).ssr.C; zeros(3,ns_ibr)] [zeros(2,ns_ibr); 1/(2*SG(3).H)*IBR(2).ssr.C; zeros(3,ns_ibr)] [ zeros(4,ns_ibr);  1/(2*SG(5).H)*IBR(3).ssr.C; zeros(1,ns_ibr)]];

Grid.ssr.A = [KD_m                          L_m                 IBR_m;
              eye(nbus)                     zeros(nbus, nbus)   zeros(nbus, nibrs*ns_ibr);
              zeros(nibrs*ns_ibr, 2*nbus)   blkdiag(IBR(1).ssr.A, IBR(2).ssr.A, IBR(3).ssr.A, IBR(4).ssr.A)];

SG_Bm   = [];

for i=1:nbus
    SG_Bm = blkdiag(SG_Bm, 1/(2*SG(i).H)) ;
end

IBR_Bm  = [zeros(nbus, nibrs*nu_ibr); blkdiag(IBR(1).ssr.B, IBR(2).ssr.B,  IBR(3).ssr.B, IBR(4).ssr.B )];

Grid.ssr.B = blkdiag(SG_Bm, IBR_Bm);

Grid.ssr.C = eye(size(Grid.ssr.A,1));

Grid.ssr.D = zeros(size(Grid.ssr.C,1), size(Grid.ssr.B,2));

stbe(Grid.ssr,0)

niter=30;

lambda_compensator = -0.01;

H = [eye(nbus) zeros(nbus, nbus) zeros(size(IBR_m,1), size(IBR_m,2))];

AG(1).Ru = 1; % Control weighted matrix of SG1
AG(2).Ru = 3; % Control weighted matrix of SG2
AG(3).Ru = 1; % Control weighted matrix of SG3
AG(4).Ru = 1; % Control weighted matrix of SG4
AG(5).Ru = 1; % Control weighted matrix of SG5
AG(6).Ru = 3; % Control weighted matrix of SG6
AG(7).Ru = 0.5*eye(2); % Control weighted matrix of IBR1
AG(8).Ru = 0.5*eye(2); % Control weighted matrix of IBR2
AG(9).Ru = 0.5*eye(2); % Control weighted matrix of IBR3
AG(10).Ru = 0.3*eye(2); % Control weighted matrix of IBR3

%Augmented parameters of the objective function for deviation system
AG(1).Q = blkdiag(eye(size(Grid.ssr.A)),1.5e4, 1e4, 0.5e4, 1.5e4, 2e4, 1e4);
AG(2).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e3, 1e3, 0.5e3, 1.5e3, 1e3, 1e3);
AG(3).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.5e4, 1.5e4, 1e4, 2e4);
AG(4).Q = blkdiag(eye(size(Grid.ssr.A)),1.5e4, 1e4, 0.5e4, 1.5e4, 1e4, 2e4);
AG(5).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.5e4, 2.5e4, 1e4, 1e4);
AG(6).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e3, 1e3, 1.5e3, 2.5e3, 1e3, 1e3);
AG(7).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.5e4, 2.5e4, 1e4, 1e4);
AG(8).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.5e4, 2.5e4, 1e4, 1e4);
AG(9).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.5e4, 2.5e4, 1e4, 1e4);
AG(10).Q = blkdiag(eye(size(Grid.ssr.A)),0.5e4, 1e4, 1.0e4, 0.5e4, 1.5e4, 1e4);

nAG = 10;
BAG = [1 1; 2 2; 3 3; 4 4; 5 5;6 6; 7 8; 9 10; 11 12; 13 14];
[pre_ssr, AG, post_ssr] = general_ncgame(Grid, AG, nAG, BAG, H, lambda_compensator, niter);

%

[P,D] = eig(post_ssr.A);
QT = inv(P);
Q = QT';
j=1;
% j is index on columns (modes)
% i is index on rows (states)
pf = [];
while j<=size(eig(post_ssr.A),1)
i=1;
while i<=size(post_ssr.A,1)
pf(i,j)=Q(i,j)*P(i,j);
i=i+1;
end
j=j+1;
end

pf_mag = round(abs(pf),2);

v = zeros(1, length(D));
vv = cell(1, length(D));
for k = 1 : length(D)
    v(k) = D(k, k);
    string_eig = string(D(k, k)); 
    vv{1,k} = string_eig(1);
end

TPF = array2table([round(v,2);pf_mag]);
%writetable(TPF, 'pf_mag_1.csv', 'Delimiter',',')



%

% t: indicates the time dimension. If the time dimension is the regular,
% then put 1. Sometimes, t_new = t * wo , then you may input wo.

function [Kp, Ki, pre_ssr, post_ssr] = LQR_Int_v1(A, B, H, t, Q, R)

ns = size(A,1); % number of states of A
nu = size(B,2); % number of inputs of B
nr = size(H,1); % number of references to track

if t == 0
    fprintf("Error: t cannot be zero, put 1 \n")
    return
end

if size(H,2) ~= ns
    fprintf("Error: The number of columns of H must be equal to the number of states of A \n")
    return
end

% Build the augmented A matrix. Note that, for a rescaled time, we use 1/t;
Ag = [A          zeros(ns,nr); 
      -H*1/t 1/t*zeros(nr,nr)];

Bg = [B ; 1/t*zeros(nr,nu)];

pre_ssr.A = Ag;
pre_ssr.B = Bg;

K = lqr(Ag, Bg, Q, R);

Kp  = K(:,1:ns);
Ki  = K(:,(ns+1):end);

A_new = (Ag - Bg * K);
B_new = [zeros(ns, nr); 1/t*eye(nr)];

post_ssr.A = A_new;
post_ssr.B = B_new;

stbe(post_ssr,0);

return
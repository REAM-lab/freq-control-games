function [est, L] = LTR_estimator(ssr, t, N)

A = ssr.A;
B = ssr.B;
C = ssr.C;

n = length(A);
Mo = eye(n);
G = eye(n);

W = t^2*(Mo) + B*B'; % Covariance of process noise
V = t^2*N ; % Covariance of measurement noise

[cov,~,~] = icare(A', ...
                         [], ...
                         (G)*W*(G)', ...
                         [], ...
                         [], ...
                         [], ...
                         -C'*V^(-1)*C); % A*P1 + P1*A' + G*W*G' - P*C'*V^(-1)*C*P

L = (cov)*(C)'*(V)^(-1);

est.A = A - (L)*(C); 
est.B = [B L]; 
est.C = eye(n,n);
est.D = zeros(n,size(est.B,2));

end
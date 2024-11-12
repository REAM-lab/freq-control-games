function [eig_vec, ssr_uc] = stbe(ssr, sprCtrb)

eig_vec = eig(ssr.A);
n = size(ssr.A, 1);

try
   q = rank(ctrb(ssr.A, ssr.B));
catch 
   fprintf('No possible to compute rank(ctrb(ssr.A, ssr.B). Zero will be shown.')
   q = 0 ;
end

cA = cond(ssr.A);
eigs_real = abs(real(eig_vec));

fprintf(' \n')
fprintf(' %19s %4.2f\n', 'Cond(A)', cA)
fprintf(' %19s %4.2f \n', 'Rank(Qctrb)', q)
fprintf(' %19s %4.2f \n', 'Max(|Re(lambda)|)', max(eigs_real))
fprintf(' %19s %4.2f \n', 'Min(|Re(lambda)|)', min(eigs_real))
fprintf(' %19s %4.2f \n', 'Ratio(Max/Min)', max(eigs_real)/min(eigs_real))

fprintf(' \n')
fprintf('  %18s  %7s  %9s  %20s  %1s \n', 'eigenvalue(A)', 'damping', 'frequency', 'rank([lambda*I-A,B])', 'n');
fprintf('%48s \n', '-----------------------------------------------------------------------------------');
for j=1:length(eig_vec)
    lambda = eig_vec(j);
    r = rank([(lambda*eye(n)-ssr.A), ssr.B]);      
    fprintf('%9.2f+%9.2fj%9.2f%11.2f%22d%3d \n', real(lambda), imag(lambda), cos(atan(imag(lambda)/real(lambda))), abs(lambda),r, n);
end

if sprCtrb==1
    [ssr_uc.Abar, ssr_uc.Bbar, ssr_uc.Cbar, ~] = ctrbf(ssr.A,ssr.B,ssr.C);
    ssr_uc.Auc = ssr_uc.Abar(1:(n-q), 1:(n-q));
    ssr_uc.Ac = ssr_uc.Abar((n-q+1):end, (n-q+1):end);
else
    ssr_uc = 0;
end

end
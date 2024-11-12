function [pre_ssr, AG, post_ssr] = general_ncgame(Grid, AG, nAG, BAG, H, lambda_compensator, niter)

ns = size(Grid.ssr.A, 1); % number of states of the grid
nr = size(H,1); % number of references to track

%Augmented state-space representation to track input
A       = [Grid.ssr.A  zeros(ns,nr); 
            -H           lambda_compensator*eye(nr,nr)];

B = [];
for i=1:nAG
    c1 = BAG(i,1);
    c2 = BAG(i,2);
    AG(i).B = [Grid.ssr.B(:,c1:c2); zeros(nr,c2-c1+1)];
    AG(i).S = (AG(i).B)*(AG(i).Ru^(-1))*(AG(i).B)';
    B= [B AG(i).B];
end

pre_ssr.A = A;
pre_ssr.B = B;

fprintf(' \n');
fprintf('[A_aug, B_aug] analysis: \n');
stbe(pre_ssr, 0);
fprintf(' \n');

for i=1:nAG
    AG(i).Ak = A;
end

fprintf('Lyapunov iterations: \n');
fprintf('  %9s  %3s  %3s \n', 'Iteration', 'Max(Ric)', 'Min(Ric)');

for k=1:niter

    for i=1:nAG
        AG(i).P = icare(AG(i).Ak,[],AG(i).Q,[],[],[],-AG(i).S);
    end

    for i=1:nAG
        sum = 0;
        for j=1:nAG
                if j~=i
                    sum = A  - (AG(j).S)*(AG(j).P);
                end
        end
        AG(i).Ak = sum ;

        AG(i).CReq = AG(i).Ak'*AG(i).P + AG(i).P*AG(i).Ak - AG(i).P*AG(i).S*AG(i).P+ AG(i).Q ;

        AG(i).F = -(AG(i).Ru)^(-1)*AG(i).B'*AG(i).P;

        % Controllers
        AG(i).Kp = AG(i).F(:,1:ns);
        AG(i).Ki = AG(i).F(:,(ns+1):end);

    end

    maxC = zeros(nAG,1);
    minC = zeros(nAG,1);
    fA = 0;

    for i=1:nAG
        maxC(i) = max(max(AG(i).CReq));
        minC(i) = min(min(AG(i).CReq));
        fA = fA + AG(i).B*AG(i).F;
    end

        fprintf('%11d %2.2e %2.2e \n', k, max(maxC), min(minC));

end

post_ssr.A = A + fA;
post_ssr.B = [];

fprintf(' \n');
fprintf('[A_aug + sum(Bi*Fi), [] ] analysis: \n');
stbe(post_ssr, 0);
fprintf(' \n');

if lambda_compensator~=0
	A_lambda_zero = [Grid.ssr.A  zeros(ns,nr); 
            		  	  -H          zeros(nr,nr)];
        fA = 0;
        for i=1:nAG
        fA = fA + AG(i).B*AG(i).F;
        end
        
        post_ssr.A_lambda_zero = A_lambda_zero + fA;
        post_ssr.B_lambda_zero = [];
 end


return

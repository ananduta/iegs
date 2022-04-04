function [p,o] = centralized_SCP_SOCP(p)


% Initialize flow solution by solving centralized SOCP with penalty
% approach
p.Gamma_pen_flag = 1;
[~,o] = solveCentralized_SOCP_he18_MA(p);

for i=1:p.n
    for jj=1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        scp.phi_old{i,j} = o.phi{i,j};
    end
end

p.Gamma_pen_flag = 0;

% Initialization of other parameters:
scp.pen_fac = 0.1;
scp.pen_fac_max = 1e5;
scp.eps_z = 1e-5;
scp.eps_s = 1e-5;
scp.v = 3;
scp.z(1) = 2;
k=1;
while 1
    % Obtain additional constraint
    p.fixedInt_flag = 0;
    [~,o,scp] = solveCentralized_SOCP_he18_MA_SCP(p,scp,k);
    
    
    er = [norm(scp.gapZ(:,k),inf);norm(scp.gapS(:,k),inf)]
    
    scp.pen_fac = min([scp.v*scp.pen_fac,scp.pen_fac_max]);
    
    if norm(scp.gapZ(:,k),inf) <= scp.eps_z && norm(scp.gapS(:,k),inf) <= scp.eps_s
        break
    end
    
    k=k+1;
end


end
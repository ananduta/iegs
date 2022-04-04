function [p,o] = centralized_SCP_SOCP_SA(p)


% Initialize flow solution by solving centralized SOCP with penalty
% approach
p.Gamma_pen_flag == 1;
[p,o] = solveCentralized_SOCP_he18(p);

for i=1:p.n
    for jj=1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        scp.phi_old{i,j} = o.phi{i,j};
    end
end

% Initialization of other parameters:
scp.pen_fac = 0.1*ones(p.nA,1);
scp.pen_fac_max = 1e3;
scp.eps_z = 1e-3;
scp.eps_s = 1e-3;
scp.v = 3;
scp.z(:,1) = 0*ones(p.nA,1);
scp.gapZ(:,1) = ones(p.nA,1);
scp.gapS(:,1) = ones(p.nA,1);
k=1;

while 1
    % Generate matrices for linear constraints       
    p = build_mat_iegs_he18(p);
    for a = 1:p.nA
        if scp.gapZ(a,k) <= scp.eps_z && scp.gapS(a,k) <= scp.eps_s
            scp.gapZ(a,k+1) = scp.gapZ(a,k);
            scp.gapS(a,k+1) = scp.gapS(a,k);
            continue
        end
        
        [p,o,scp] = solveCentralized_SOCP_he18_SA_SCP(p,scp,a);

        

        scp.pen_fac(a) = min([scp.v*scp.pen_fac(a),scp.pen_fac_max]);
    end
    
    k=k+1;
    if norm(scp.gapZ(:,k),inf)<= scp.eps_z && norm(scp.gapS(:,k),inf)<= scp.eps_S
        break
    end
end


end
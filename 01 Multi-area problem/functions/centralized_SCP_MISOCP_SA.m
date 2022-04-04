function [p,o] = centralized_SCP_MISOCP_SA(p)

% Initialization
%% Identify dimensions of decision variables (per time step, h)
    
    for i = 1:p.n
        p.nxc(i) = length(p.en.Nex{i});   % auxiliary consensus variables
        p.nx(i) = 5+p.en.noN(i)+p.nxc(i);
        p.nt(i) = p.gn.noN(i);
        p.ny(i) = 2+2*p.gn.noN(i)+p.nt(i);
        p.nz(i) = p.gn.noN(i);
        p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
        
        p.id_phi{i} = 0;
        for jj = 1:p.gn.noN(i)
            p.id_phi{i}(jj) = p.nx(i) + 2+(2)*(jj-1)+1;
        end
    end
    p.id_dg = 1;
    p.id_ch = 2;
    p.id_dh = 3;
    p.id_gu = 4;
    p.id_th = 5;
    % y vars. must +p.nx(i)
    p.id_psi = 1;
    p.id_gs = 2;

% Generate matrices for linear constraints       
p = build_mat_iegs_he18(p);
    
% Initialize flow solution by solving centralized MISOCP with penalty
% approach

for a = 1:p.nA
    p.Gamma_pen_flag = 1;
    [~,o] = solveCentralized_MISOCP_he18_SA(p,a);
    for ii=1:length(p.Ae{a})
        i = p.Ae{a}(ii);
        for jj=1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
            scp.phi_old{i,j} = o.phi{i,j};
        end
    end
end

p.Gamma_pen_flag = 0;


% Initialization of other parameters:
scp.pen_fac = 0.1*ones(p.nA,1);
scp.pen_fac_max = 1e6;
scp.eps_z = 1e-3;
scp.eps_s = 1e-3;
scp.v = 3;
scp.z(:,1) = 0*ones(p.nA,1);
scp.gapZ(:,1) = ones(p.nA,1);
scp.gapS(:,1) = ones(p.nA,1);
scp.gapZ(:,2) = ones(p.nA,1);
scp.gapS(:,2) = ones(p.nA,1);
k=1;

while 1
    
    
    
    
    for a = 1:p.nA
        if scp.gapZ(a,k) <= scp.eps_z && scp.gapS(a,k) <= scp.eps_s
            scp.gapZ(a,k+1) = scp.gapZ(a,k);
            scp.gapS(a,k+1) = scp.gapS(a,k);
            continue
        end
        
        [~,oo,scp] = solveCentralized_MISOCP_he18_SA_SCP(p,scp,a,k);
        for ii = 1:length(p.Ae{a})
            i = p.Ae{a}(ii);
            o.u{i} = oo.u{i};

            o.p_di{i} = p.m.Sdg{i}*o.u{i};
            o.p_ch{i} = p.m.Sch{i}*o.u{i};
            o.p_dh{i} = p.m.Sdh{i}*o.u{i};
            %o.p_mg{i} = p.m.Smg{i}*o.u{i};
            o.d_gu{i} = p.m.Sgu{i}*o.u{i};
            o.th{i} = p.m.Sth{i}*o.u{i}/p.scaling;
            %o.v{i} = p.m.Sv{i}*o.u{i}/p.scaling;
            for jj = 1:length(p.en.N{i})
                j = p.en.N{i}(jj);
                o.p_l{i,j} = p.m.Spl{i,j}*o.u{i};
            end

        % y
            o.psi{i} = p.m.Spsi{i}*o.u{i};
            o.gs{i} = p.m.Sgs{i}*o.u{i};


            for jj = 1:length(p.gn.Nin{i})
                j = p.gn.Nin{i}(jj);
                o.phi{i,j} = p.m.Sphi{i,j}*o.u{i};
                o.pen{i,j} = p.m.Spen{i,j}*o.u{i};

                o.alpha{i,j} = p.m.Salp{i,j}*o.u{i};
            end
            for jj = 1:length(p.gn.Nex{i})
                j = p.gn.Nex{i}(jj);
                o.phi{i,j} = p.m.Sphi{i,j}*o.u{i};
                o.pen{i,j} = p.m.Spen{i,j}*o.u{i};


            end
        end
        

        scp.pen_fac(a) = min([scp.v*scp.pen_fac(a),scp.pen_fac_max]);
    end
    
    er = [norm(scp.gapZ(:,k),inf);norm(scp.gapS(:,k),inf)]
    if norm(scp.gapZ(:,k),inf)<= scp.eps_z && norm(scp.gapS(:,k),inf)<= scp.eps_s
        break
    end
    k=k+1;
end


end
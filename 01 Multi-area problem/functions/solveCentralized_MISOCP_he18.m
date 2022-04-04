function [p,o]= solveCentralized_MISOCP_he18(p)
% Computing v-GNE of Mixed-integer IEGS
% W. Ananduta
% 17/02/2022

    
    %% Identify dimensions of decision variables (per time step, h)
    for i = 1:p.n
        p.nx(i) = 5+p.en.noN(i);
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


    % Generate agents' indices in the concatenated decision variable u
    p.idAgent = index_decision(p.nu,p.n,p.h);
    idStart = p.idAgent(:,1);
    idEnd = p.idAgent(:,2);
       
    % Generate matrices for linear constraints       
    p = build_mat_iegs_he18(p);  
    
    % local constraints
    A_all = [];
    b_all = [];
    for i = 1:p.n
        A = p.m.Aineq{i};
        b = p.m.bineq{i};
        Aeq = p.m.Aeq{i};
        beq = p.m.beq{i};
        
        At{i} = [A;Aeq;-Aeq];
        bt{i} = [b;beq;-beq];
        
        Ai = zeros(size(At{i},1),idEnd(p.n));
        Ai(:,idStart(i):idEnd(i)) = At{i};
        bi = bt{i};
        
        A_all = sparse([A_all; Ai]);
        b_all = [b_all; bi];
    end
    
    % coupling constraints
    for i = 1:p.n
        for jj=1:length(p.en.N{i})
            j = p.en.N{i}(jj);
            
            %power flow constraints
            Ac = zeros(p.h,idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = p.m.PF{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.PF{i,j}{2};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
        end
        
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
                        
            % gas-flow coupling constraints
            Ac = zeros(size(p.m.Gc{i,j}{1},1),idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = p.m.Gc{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.Gc{i,j}{2};
            
            bc = p.m.gc{i,j};
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
            
        end
        
                
    end
        
    model.A = A_all;
    model.rhs = b_all;
    
    % Quadratic constraints from coupling between d_gu and p_gu
    countQC = 1;
    for i = 1:p.n
        if p.en.dgu_un(i) == 1
            for h = 1:p.h
                id_p_dg_ih = idStart(i)-1+p.nu(i)*(h-1)+p.id_dg;
                id_p_gu_ih = idStart(i)-1+p.nu(i)*(h-1)+p.id_gu;
                
                Qc = sparse(zeros(idEnd(p.n)));
                Qc(id_p_dg_ih,id_p_dg_ih) = p.en.q1(i);
                
                q = sparse(zeros(idEnd(p.n),1));
                q(id_p_dg_ih,1) = p.en.q2(i);
                q(id_p_gu_ih,1) = -1;
                
                rhs = -p.en.q3(i);
                
                model.quadcon(countQC).Qc = Qc;
                model.quadcon(countQC).q = q;
                model.quadcon(countQC).rhs = rhs;
                countQC = countQC+1;
            end
        end
    end
    
    % Quadratic constraints from gas-flow approximation
    
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
            
            for h = 1:p.h
                id_gamma =  p.nx(i) + 2+(2)*jj;
                id_gamma_ijh = idStart(i)-1+p.nu(i)*(h-1)+p.nx(i) + 2+(2)*jj;
                id_phi_ijh = idStart(i)-1+p.nu(i)*(h-1)+ p.id_phi{i}(jj);
                
                Qc = sparse(zeros(idEnd(p.n)));
                Qc(id_phi_ijh,id_phi_ijh) = 1/p.gn.cf(i,j);
                
                q = sparse(zeros(idEnd(p.n),1));
                q(id_gamma_ijh) = -1;
                
                rhs = 0;
                
                model.quadcon(countQC).Qc = Qc;
                model.quadcon(countQC).q = q;
                model.quadcon(countQC).rhs = rhs;
                countQC = countQC+1;
            end
        end
    end
    
    % Generate matrices for cost function
    
    p = alg_param5(p); % nothing is used from the output of this function, but must be run.
    
    p = build_mat_cost_pen(p);
    
    % linear term
    f = cat(1,p.m.ch{:});
    
    % quadratic term
    
    Hloc = sparse(zeros(idEnd(p.n)));
    
    for i = 1:p.n
        Hloc(idStart(i):idEnd(i),idStart(i):idEnd(i)) = 2*p.m.Qph{i};
    end
    
    
   % p.m.Hloc = Hloc;
    
   % p.m.Hall = sparse(p.m.Hloc); 
    
   % H = p.m.Hall; 
    
    %J = x'*Hloc*x + f'*x;
    model.Q = Hloc;
    model.obj = f;
    
    
    % lower bound
    model.lb = -inf*ones(idEnd(p.n),1);
    
    
    % set continuous and binary variables
    for i = 1:p.n
%         for j = idStart(i):idEnd(i)
%             if j <= idStart(i) + p.nx(i)+p.ny(i)
%                 model.vtype(j,:) = 'C';
%             else
%                 model.vtype(j,:) = 'B';
%             end
%             
%         end
        for j =1:p.nu(i)
            if j <=  p.nx(i)+p.ny(i)
                vtype{i}(j,:) = 'C';
            else
                vtype{i}(j,:) = 'B';
            end
        end
    end
    model.vtype = cat(1,vtype{:});
    %% solve with Gurobi
    

    res = gurobi(model);
    o.res = res;
    u = res.x;
    o.u_all= u;
%    o.res = res;
    
    
    % Assign solution
    for i = 1:p.n
        o.u{i} = u(idStart(i):idEnd(i),1);

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
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            o.phi{i,j} = p.m.Sphi{i,j}*o.u{i};
            o.pen{i,j} = p.m.Spen{i,j}*o.u{i};
        end
    end
    
end
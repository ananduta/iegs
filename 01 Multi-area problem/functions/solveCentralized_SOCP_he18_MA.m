function [p,o]= solveCentralized_SOCP_he18_MA(p)
% Computing v-GNE of convexified game centrally via minimizing potential function
% W. Ananduta
% 07/02/2022
% Updated for Convex problem 15/02/2022

    
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
        
        Ai = sparse(zeros(size(At{i},1),idEnd(p.n)));
        Ai(:,idStart(i):idEnd(i)) = At{i};
        bi = bt{i};
        
        A_all = sparse([A_all; Ai]);
        b_all = [b_all; bi];
    end
    
    % coupling constraints
    for i = 1:p.n
        for jj=1:length(p.en.Nin{i})
            j = p.en.Nin{i}(jj);
            
            %power flow constraints (local)
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.PF{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.PF{i,j}{2};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
        end
        
        for jj=1:length(p.en.Nex{i})
            j = p.en.Nex{i}(jj);
            
            %if ~isempty(find(p.en.Nex{i}==j))
                %power flow constraints (local) with consensus aux. var
                % power flow equation


                PFi = zeros(1,p.nu(i));
                PFi(1,p.id_th) = p.en.Bnet(i,j);
                %PFi(1,p.id_v) = -p.en.Gnet(i,j);
                PFi(1,p.id_th+jj) = -1;
                PFi(1,p.nx(i)-p.nxc(i)+jj) = -p.en.Bnet(i,j);
                PFit = kron(eye(p.h),PFi);

                Ac = sparse(zeros(p.h,idEnd(p.n)));
                Ac(:,idStart(i):idEnd(i)) = PFit;
                bc = zeros(p.h,1);

                A_all = sparse([A_all;Ac]);
                b_all = [b_all;bc];


                %consensus constraints on theta (with other Areas)
                Ac = sparse(zeros(p.h,idEnd(p.n)));
                Ac(:,idStart(i):idEnd(i)) = p.m.Sthc{i,j};
                Ac(:,idStart(j):idEnd(j)) = -p.m.Sth{j};

                Ac = [Ac;-Ac];
                bc = zeros(2*p.h,1);

                A_all = sparse([A_all;Ac]);
                b_all = [b_all;bc];
            %end
        end
        
        for jj = 1:length(p.gn.Nex{i})
            j = p.gn.Nex{i}(jj);
            
            % reciprocity constraints
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.Sphi{i,j};
            Ac(:,idStart(j):idEnd(j)) = p.m.Sphi{j,i};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
        end
        
        for jj = 1:length(p.gn.Nin{i})
            j = p.gn.Nin{i}(jj);
            
            
             % reciprocity constraints
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.Sphi{i,j};
            Ac(:,idStart(j):idEnd(j)) = p.m.Sphi{j,i};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
%             % agreement on binary vars
%             Ac = sparse(zeros(p.h,idEnd(end)));
%             Ac(:,idStart(i):idEnd(i)) = p.m.Salp{i,j};
%             Ac(:,idStart(j):idEnd(j)) = p.m.Salp{j,i};
%             bc = ones(p.h,1);
%             
%             Ac = [Ac;-Ac];
%             bc = [bc;-bc];
%             
%             A_all = sparse([A_all;Ac]);
%             b_all = [b_all;bc];
            
            % gas-flow coupling constraints
            Ac = sparse(zeros(size(p.m.Gc{i,j}{1},1),idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.Gc{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.Gc{i,j}{2};
            
            bc = p.m.gc{i,j};
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
            % Additional constraints if binary variables are fixed
            if p.fixedInt_flag == 1
                Ac = sparse(zeros(p.h,idEnd(p.n)));
                Ac(:,idStart(i):idEnd(i)) = p.m.Salp{i,j};
                
                bc = p.alphaFixed{i,j};
                
                Ac = [Ac;-Ac];
                bc = [bc;-bc];
                
                A_all = sparse([A_all;Ac]);
                b_all = [b_all;bc];
                
            end
            
            
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
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            if ~isempty(find(p.gn.Nin{i}==j))
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
    %% solve with Gurobi
    
    params.outputflag = 1;
    params.BarHomogeneous = 1;
    params.DualReductions = 0;
    params.NumericFocus = 3;
    res = gurobi(model,params);
    if res.status(1) == 'N' || res.status(1) == 'I'
        o.flag_NumIssue = 1;
    else
        o.res = res;
        u = res.x;
        o.u_all= u;
        o.flag_NumIssue = 0;
    %    o.res = res;


        % Assign solution
        for i = 1:p.n
            o.u{i} = u(idStart(i):idEnd(i),1);

            o.p_di{i} = p.m.Sdg{i}*o.u{i};
            o.p_ch{i} = p.m.Sch{i}*o.u{i};
            o.p_dh{i} = p.m.Sdh{i}*o.u{i};
            %o.p_mg{i} = p.m.Smg{i}*o.u{i};
            o.d_gu{i} = p.m.Sgu{i}*o.u{i};
            o.th{i} = p.m.Sth{i}*o.u{i};
            %o.v{i} = p.m.Sv{i}*o.u{i}/p.scaling;
            for jj = 1:length(p.en.Nin{i})
                j = p.en.Nin{i}(jj);
                o.p_l{i,j} = p.m.Spl{i,j}*o.u{i};

            end

            for jj = 1:length(p.en.Nex{i})
                j = p.en.Nex{i}(jj);
                o.p_l{i,j} = p.m.Spl{i,j}*o.u{i};

                o.thc{i,j} = p.m.Sthc{i,j}*o.u{i};


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
    end
    
end
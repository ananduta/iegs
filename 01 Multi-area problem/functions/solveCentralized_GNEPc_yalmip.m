function [p,o]= solveCentralized_GNEPc_yalmip(p)
% Computing v-GNE of convexified game centrally via minimizing potential function
% W. Ananduta
% 07/02/2022
% Updated for Convex problem 15/02/2022

    
    % Generate agents' indices in the concatenated decision variable u
    p.idAgent = index_decision(p.nu,p.n,p.h);
    idStart = p.idAgent(:,1);
    idEnd = p.idAgent(:,2);
    
    % Initializing yalmip environment
    x = sdpvar(idEnd(end),1);
    J = 0;
    cons = [];
    
    % Generate matrices for  constraints       
    p = build_mat_iegs_opt(p);  
    
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
        b_all = sparse([b_all; bi]);
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
            b_all = sparse([b_all;bc]);
        end
        
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
            % reciprocity constraints
            Ac = zeros(p.h,idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = p.m.Sphi{i,j};
            Ac(:,idStart(j):idEnd(j)) = p.m.Sphi{j,i};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = sparse([b_all;bc]);
            
            % gas flow constraints
            Ac = zeros(p.h,idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = p.m.Hc{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.Hc{i,j}{2};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = sparse([b_all;bc]);
            
            % psi coupling constraints
            Ac = zeros(size(p.m.G{i,j}{1},1),idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = p.m.G{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.G{i,j}{2};
            
            bc = p.m.g{i,j};
            
            A_all = sparse([A_all;Ac]);
            b_all = sparse([b_all;bc]);
            
            
        end
        
                
    end
        
%    p.m.A_all = A_all;
%    p.m.b_all = b_all;
    
    
    cons = [cons, A_all*x <= b_all];
    
    % Quadratic constraints from coupling between d_gu and p_gu
    for i = 1:p.n
        if p.en.dgu_un(i) == 1
            for h = 1:p.h
                p_dg_h = x(idStart(i)-1+p.nu(i)*(h-1)+p.id_dg);
                p_gu_h = x(idStart(i)-1+p.nu(i)*(h-1)+p.id_gu);
                cons = [cons, p.en.q1(i)*p_dg_h^2 +p.en.q2(i)*p_dg_h + p.en.q3(i) - p_gu_h <=0];
            end
        end
    end
    
    
    % Generate matrices for cost function
    
    p = alg_param4(p); % nothing is used from the output of this function, but must be run.
    
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
    
    J = x'*Hloc*x + f'*x;
    
    %% solve with Yalmip
    

    
    options = sdpsettings('verbose',1,'solver','fmincon');
    %options = sdpsettings('verbose',1,'solver','bnb');
    sol = optimize(cons,J,options);
    
    u = value(x);
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
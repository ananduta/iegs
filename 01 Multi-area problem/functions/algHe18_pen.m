function [p,o,q] = algHe18_pen(p)
% Algorithm 3 of He et al TPS'18 Decentralized...
% W. Ananduta
% 14/02/2022


%% Step 3.1 INITIALIZATION
% Solve convexified SOCP multi-area problem
 p.fixedInt_flag =0;
% CENTRALIZED
% with SCP

% [~,o] = centralized_SCP_SOCP(p);


% % with penalty
tic
p.Gamma_pen_flag = 1;
[~,o] = solveCentralized_SOCP_he18_MA(p);
q.time(1) = toc;
% (OR) DECENTRALIZED via ADMM
% with SCP

% with penalty

BinVarsOld = cat(1,o.alpha{:});
% Gather initial value of coupling variables/shared information

%% ITERATION
rmax = 30;
for r = 1:rmax
    %% Step 3.2
    tic
    % Averaging shared information/value of coupling variables
    for i = 1:p.n
        p.thc{i} = [];
        if ~isempty(p.en.Nex{i})
            p.thc{i} = o.th{i};
            for jj = 1:length(p.en.Nex{i})
                j = p.en.Nex{i}(jj);
                p.thc{i} =  p.thc{i} + o.thc{j,i};
            end
            p.thc{i} = p.thc{i}/(length(p.en.Nex{i})+1);
        end

        for jj=1:length(p.gn.Nex{i})
            j = p.gn.Nex{i}(jj);
            p.phic{i,j} = 0.5*(o.phi{i,j}-o.phi{j,i});
        end

    end

    % Solve single area MI-SOCP problem

    % with SCP
%     [~,o] = centralized_SCP_MISOCP_SA(p);


    % with penalty-based approach
    [~,o] = centralized_pen_SA(p);

    p.alphaFixed = o.alpha;
    BinVars = cat(1,o.alpha{:});
    %% Step 3.3
    % Solve convex SOCP multi-area problem with fixed integer
    p.fixedInt_flag =1;
    
    % CENTRALIZED
    % with SCP
    
%     [~,o] = centralized_SCP_SOCP(p);
    
    % with penalty
    p.Gamma_pen_flag = 1;
    [p,o]= solveCentralized_SOCP_he18_MA(p);
    q.time(r+1) = toc;
    if o.flag_NumIssue == 1
        q = 0;
        break
    end
    % (OR) DECENTRALIZED via ADMM
    % with SCP

    % with penalty


    %% Step 3.4
    % Evaluate cost
    for i=1:p.n
    [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
    end
    q.Jt(r) = sum(q.J(:,r));
    q.Pt(r) = sum(q.P(:,r));
    q.er_gf = gasFlow_error(p,o);
    % Check stopping criterion
    q.gapI(:,r) = BinVars - BinVarsOld;
    BinVarsOld = BinVars;
    if norm(q.gapI(:,r),inf) <= 1e-7
        break
    end
    %r=r+1;
end


end
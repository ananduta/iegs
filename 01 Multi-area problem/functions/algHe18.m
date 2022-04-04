function [p,o] = algHe18(p)
% Algorithm 3 of He et al TPS'18 Decentralized...
% W. Ananduta
% 14/02/2022


%% Step 3.1 INITIALIZATION
% Solve convexified SOCP multi-area problem
 p.fixedInt_flag =0;
% CENTRALIZED
% with SCP

[~,o] = centralized_SCP_SOCP(p);


% % with penalty
% p.Gamma_pen_flag = 1;
% [~,o] = solveCentralized_SOCP_he18_MA(p);

% (OR) DECENTRALIZED via ADMM
% with SCP

% with penalty

BinVarsOld = cat(1,o.alpha{:});
% Gather initial value of coupling variables/shared information

%% ITERATION
while 1
    %% Step 3.2

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
    [~,o] = centralized_SCP_MISOCP_SA(p);


    % with penalty-based approach
%     [~,o] = centralized_pen_SA(p);
    
    p.alphaFixed = o.alpha;
    BinVars = cat(1,o.alpha{:});
    %% Step 3.3
    % Solve convex SOCP multi-area problem with fixed integer
    p.fixedInt_flag =1;
    
    % CENTRALIZED
    % with SCP
    
    [~,o] = centralized_SCP_SOCP(p);
    
    % with penalty
%     p.Gamma_pen_flag = 1;
%     [~,o]= solveCentralized_SOCP_he18_MA(p);


    % (OR) DECENTRALIZED via ADMM
    % with SCP

    % with penalty


    %% Step 3.4
    % Check stopping criterion
    gapI = BinVars - BinVarsOld;
    BinVarsOld = BinVars;
    if norm(gapI,inf) == 0
        break
    end
end

end
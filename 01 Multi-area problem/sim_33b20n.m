% Simulations
% P2P market on IEGDS
% 05/08/2021


clear all
close all
clc

run('pathdef')

% Add path of folder 'functions'
addpath([pwd,'/functions'])


% generate case
ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 33;
%n_passive = 0;


run('gen_iegs_33b_20n.m')




gamma(1) = 0;
gammaUpper = inf;
gammaLower = gamma(1);


r_max = 5;
%%
for r = 1:r_max
    %clearvars('o','e')
    
    p.pen = gamma(r);

    %% Stage 1
    % Solve convexified GNEP via proximal-point method
    [p,s] = initializeSolve_GNEPc(p,r);
    [s,p,o] = iterate_GNEPc_pen(s,p);
%     [p,o]= solveCentralized_GNEPc(p);
    save(['sim_33b20n_0702_stage1_R'],'o','p')
    %% Stage 2
    o = solve_binary(o,p);
    
    [o,e] = solve_pressure3(o,p);
    
    o.pen_w(r) = gamma(r);
    er(r) = o.gfv_max;
%    n_iter(r) = o.n_iter;
    
    
    %% Update gamma
    [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);
    
    % compute cost
    for i=1:p.n
        [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
    end
    q.Jt(r) = sum(q.J(:,r));
    q.Pt(r) = sum(q.P(:,r));

    q.gamma = gamma;
    q.er = er;
    %q.nIter = n_iter;
    q.gammaLower = gammaLower;
    q.gammaUpper = gammaUpper;
    save(['sim_33b20n_0702_R'],'o','e','p','q')
    save(['sim_33b20n_0702_R_r',num2str(r)],'o','p')
    %%
    if er(1) < 1e-5
        break
    end
end

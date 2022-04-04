% Simulations
% P2P market on IEGDS
% 05/08/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'/functions'])

%load('sim0911.mat')

% generate case
ty = [1]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 6;
n_passive = 0;

sp_set = 0.7;
run('gen_iegs_6n.m')




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

    %% Stage 2
    o = solve_binary(o,p);
    
    [o,e] = solve_pressure3(o,p);
    
    o.pen_w(r) = gamma(r);
    er(r) = o.gfv_max;
    n_iter(r) = o.n_iter;
    if er(1) == 0
        break
    end
    %% Update gamma
    [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);
    
    % compute cost
    for i=1:p.n
        q.J(i,r) = cost_compute(o,p,i);
    end
    
    
end
%%
q.Jt = sum(q.J);

q.gamma = gamma;
q.er = er;
q.nIter = n_iter;
q.gammaLower = gammaLower;
q.gammaUpper = gammaUpper;

figure
subplot(4,1,1)
plot(q.Jt/q.Jt(1),'o-','LineWidth',1.4)
title('Total cost')
subplot(4,1,2)
plot(q.er,'o-','LineWidth',1.4)
title('Max error of gas flow equation')
subplot(4,1,3)
plot(q.gamma,'o-','LineWidth',1.4)
title('Penalty weight')
subplot(4,1,4)
plot(q.nIter,'o-','LineWidth',1.4)
title('Number of iterations')
save(['sim_6nodes_0202b'],'o','e','p','q')

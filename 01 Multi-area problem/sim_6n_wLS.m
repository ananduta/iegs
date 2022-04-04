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
ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 6;
n_passive = 0;

sp_set = 0.7;
run('gen_iegs_6n.m')




gamma(1) = 2e4;
gammaUpper = inf;
gammaLower = gamma(1);


r_max = 5;
%%
for r = 1:r_max
    %clearvars('o','e')
    
    p.pen = gamma(r);

    %% Stage 1
    % Solve convexified GNEP via proximal-point method
    [s,p,o] = solve_GNEPc_pen(p);

    %% Stage 2
    o = solve_binary(o,p);
    
    [o,e] = solve_pressure2(o,p);
    
    o.pen_w(r) = gamma(r);
    er(r) = o.gfv_max;
    n_iter(r) = o.n_iter;
    %% Update gamma
    [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper);
    
end

save(['sim_1611_wLS_a'],'o','e','er','gamma')

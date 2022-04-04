% Simulations
% P2P market on IEGDS
% 05/08/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'/functions'])


% ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
% tc = [1]; %uniform trading cost
% 
% % set the number of agents
% n_agents = 6;
% n_passive = 0;
% 
% sp_set = 0.7;
% 
% % generate case
% run('gen_iegs_6n.m')
% 
% [p.phi_max_a,p.flag] = restrict_flow(p.gn);
% % % solve with two-stage approach
% % [s,p] = TS_ED(p);

pen_v = [1];
for c = 1:length(pen_v)
    pen = pen_v(c);
    load('sim0911.mat')
    clearvars('o','e')
    for i = 1:p.n
    p.nx(i) = 4+ p.en.noN(i);
    p.nt(i) = p.gn.noN(i);
    p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
    p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
    p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
    end
    p.pen = pen;

    %% Stage 1
    % Solve convexified GNEP via proximal-point method
    [s,p,o] = solve_GNEPc_pen(p);

    o = solve_binary(o,p);
    
    [o,e] = solve_pressure2(o,p);
    o.pen_w = pen;
    save(['sim_0911_',num2str(pen)],'o','e')
end

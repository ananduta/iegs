% Simulations
% P2P market on IEGDS
% 05/08/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'/functions'])


ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 6;
n_passive = 0;

sp_set = 0.7;

% generate case
run('gen_iegs_6n.m')

[p.phi_max_a,p.flag] = restrict_flow(p.gn);
% % solve with two-stage approach
% [s,p] = TS_ED(p);

%% Stage 1
% Solve convexified GNEP via proximal-point method
[s,p,o] = solve_GNEPc(p);

o = solve_binary(o,p);

[o,e] = solve_pressure2(o,p);

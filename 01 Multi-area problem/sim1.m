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


    
run('gen_MA.m')

%load('case_73n_0203.mat')
%save('case_73n_0203.mat','p');

gamma_pen = [0:5:20];
pwa_reg = [20:20:100];
p_data = p;



for rr = 1:length(pwa_reg)
    pW = p_data;
    pW.r = pwa_reg(rr);
    p2{rr} = pW;
    [~,o2{rr},q2{rr}] = algWicak(pW);
end
    
for rr = 1:length(gamma_pen)
    pH = p_data;
    pH.gamma_pen = gamma_pen(rr);
    p1{rr} = pH;
    [~,o1{rr},q1{rr}] = algHe18_pen(pH);
end

save('sim_12n_0303b.mat','o1','o2','p1','p2','q1','q2')
        
    


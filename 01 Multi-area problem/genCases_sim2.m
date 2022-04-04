% generate 10 cases


% Add path of folder 'functions'
addpath([pwd,'/functions'])
clear all
    close all
    clc
    ndata = 10;
for d = 1:ndata
      

    %load('sim0911.mat')

    % generate case
    ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
    tc = [1]; %uniform trading cost

    % set the number of agents
    n_agents = 33;
    n_passive = 0;

    r_max = 30;
    sp_set = 0.7;

    %run('gen_iegs_6n.m')
    run('gen_iegs_33b_20n.m')
    q.gamma = 0;
    q.er = 0;
    save(['sim2_1402_b',num2str(d),'.mat'],'p','q')
    clearvars('p')
end
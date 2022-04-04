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
n_agents = 73;
%n_passive = 0;
%%
rmax = 100;
r =55 ;
while r <= rmax    
    run('gen_MA_ran.m')
    p_data = p;
        
    gamma_pen = [0,50];
    pwa_reg = [20:20:80];
    
    pH = p_data;
    pH.gamma_pen = gamma_pen(1);
    [~,o1{r},q1{r}] = algHe18_pen(pH);
    if o1{r}.flag_NumIssue == 1
       continue
    end
    
    pH.gamma_pen = gamma_pen(2);
    [~,o2{r},q2{r}] = algHe18_pen(pH);
    if o2{r}.flag_NumIssue == 1
       continue
    end
    for rr=1:length(pwa_reg)
        pW = p_data;
        pW.r = pwa_reg(rr);
        [~,o3{r,rr},q3{r,rr}] = algWicak(pW);
        save(['sim5_472n_1903.mat'],'q1','q2','q3')
    end


    
    r = r+1;
end
        
    


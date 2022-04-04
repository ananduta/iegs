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

rmax = 10;
r =7 ;
while r <= rmax    
    run('gen_MA_ran.m')
    p_data = p;

    

    gamma_pen = [0:20:100];
    pwa_reg = [20:20:80];

    for rr = 1:length(gamma_pen)
        pH = p_data;
        pH.gamma_pen = gamma_pen(rr);
        p1{rr} = pH;
        [~,o1{rr},q1{rr}] = algHe18_pen(pH);
        if o1{rr}.flag_NumIssue == 1
            break
        end
            
    end
    if o1{rr}.flag_NumIssue == 1
       continue
    end
    for rr = 1:length(pwa_reg)
        pW = p_data;
        pW.r = pwa_reg(rr);
        p2{rr} = pW;
        [~,o2{rr},q2{rr}] = algWicak(pW);
    end


    save(['case2_472n_0903_',num2str(r),'.mat'],'p');
    save(['sim2_472n_0903_',num2str(r),'.mat'],'o1','o2','p1','p2','q1','q2')
    r = r+1;
end
        
    


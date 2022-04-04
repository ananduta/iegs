% Simulations
% P2P market on IEGDS
% test centralization
% 07/02/2021


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

rr_max = 3;
sp_set = 0.7;
for rr = 1:rr_max
    
    run('gen_iegs_33b_20n.m')




    gamma(1) = 0;
    gammaUpper = inf;
    gammaLower = gamma(1);


    r_max = 2;
    %%
    for r = 1:r_max
        %clearvars('o','e')

        p.pen = gamma(r);

        %% Stage 1
        [p,o1]= solveCentralized_GNEPc(p);

        % Solve convexified GNEP via proximal-point method
        [p,s] = initializeSolve_GNEPc(p,r);
        [s,p,o] = iterate_GNEPc_pen(s,p);

        %% Stage 2
        o = solve_binary(o,p);

        [o,e] = solve_pressure3(o,p);

        o1 = solve_binary(o1,p);

        [o1,e] = solve_pressure3(o1,p);

        o.pen_w(r) = gamma(r);
        o1.pen_w(r) = gamma(r);

        er(r) = o.gfv_max;
        er1(r) = o1.gfv_max;

        n_iter(r) = o.n_iter;
        if er(1) == 0
            break
        end
        %% Update gamma
        [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);

        % compute cost
        for i=1:p.n
            [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
            [q1.J(i,r),q1.P(i,r)] = cost_compute(o1,p,i);
        end
        q.Jt(r) = sum(q.J(:,r));
        q.Pt(r) = sum(q.P(:,r));
        q1.Jt(r) = sum(q1.J(:,r));
        q1.Pt(r) = sum(q1.P(:,r));


        q.gamma = gamma;
        q.er = er;
        q1.er = er1;

        t.diffP(rr,r) = q.Pt(r) - q1.Pt(r);
        t.Pt_dist(rr,r) = q.Pt(r);
        t.Pt_cent(rr,r) = q1.Pt(r);
        t.er_dist(rr,r) = er(r);
        t.er_cent(rr,r) = er1(r);
        t.diffU(rr,r) = norm(q1.u_all - cat(1,o.u{:}));
        %q.nIter = n_iter;
        q.gammaLower = gammaLower;
        q.gammaUpper = gammaUpper;
        save(['sim_33n_0702_R_rr',num2str(rr)],'q','q1','t')
        save(['sim_33n_0702_R_rr',num2str(rr),'_r',num2str(r)],'o','p','o1','q','q1')
        %%
        if er(1) < 1e-5
            break
        end

    end
    
end


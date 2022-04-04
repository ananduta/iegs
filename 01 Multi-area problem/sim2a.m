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
n_agents = 33;
n_passive = 0;

r_max = 30;
sp_set = 0.7;

%run('gen_iegs_6n.m')
%run('gen_iegs_33b_20n.m')

d = 5;
dat = ['sim2_1002_b',num2str(d),'.mat'];
load(dat)


%pwaReg = [10,20,30,40,50];
pwaReg1 = [10:10:50];
pwaReg2 = [15:10:45];
pwaReg3 = [55:5:60];
pwaReg = [pwaReg1 pwaReg2 pwaReg3];
%%
r_max = 30;
for rr = 1:length(pwaReg)
    
    p.r = pwaReg(rr);
        p.gn.r = p.r;

        %% Identify dimensions of decision variables (per time step, h)
       for i = 1:p.n
            p.nx(i) = 8+ p.en.noN(i);
            p.nt(i) = p.gn.noN(i);
            p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
            p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
            p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
       end
    
    gamma = q.gamma;
    er = q.er;
    
    gamma(rr,1) = 0;
    gammaUpper = inf;
    gammaLower = gamma(rr,1);
    
    %%
    for r = 1:r_max
        p.pen = gamma(rr,r);
        

        %% Stage 1
        [p,o]= solveCentralized_GNEPc(p);

        %% Stage 2
        o = solve_binary(o,p);

        [o,e] = solve_pressure3(o,p);

        

        o.pen_w(rr,r) = p.pen;
        
        er(rr,r) = o.gfv_max;
        

        % compute cost
        for i=1:p.n
            [q.J{r}(i,rr),q.P{r}(i,rr)] = cost_compute(o,p,i);
            
        end
        q.Jt(rr,r) = sum(q.J{r}(:,rr));
        q.Pt(rr,r) = sum(q.P{r}(:,rr));
        


        q.gamma = gamma;
        q.er = er;
        

       
        
        %% Update gamma
        [gamma(rr,r+1),gammaLower,gammaUpper]=gammaRules(er(rr,r),gamma(rr,r),gammaLower,gammaUpper,r);
        %q.nIter = n_iter;
        q.gammaLower = gammaLower;
        q.gammaUpper = gammaUpper;
        save(['sim2_1102_b',num2str(d)],'q')
        %save(['sim2_0902_res',num2str(r),'_rr',num2str(rr)],'o','p','o1','q','q1')
        %%
%         if er(1) < 1e-5
%             break
%         end

    end
    
    
end


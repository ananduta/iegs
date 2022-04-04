function [p,o,q] = algWicak(p)

    % Set number of regions for PWA gas-flow model
    %p.r = 30
    p.gn.r = p.r;
    
    gamma(1) = 0;
    gammaUpper = inf;
    gammaLower = gamma(1);
    

    r_max = 1;
    
    % Default parameters for some algorithms
    p.pen = 0; % for our algorithm
    p.Gamma_pen_flag = 0; % for he'18, if penalty-based used, then set to 1
    tic
    p = initialize_GNEPc_MA(p);
    q.time(1) = toc;
    %%
    for r = 1:r_max
        %clearvars('o','e')

        p.pen = gamma(r);

        %% Stage 1
        tic
    %     % Solve convexified GNEP via proximal-point method
    %     [p,s] = initializeSolve_GNEPc(p,r);
    %     [s,p,o] = iterate_GNEPc_pen(s,p);
        [p,o]= solveCentralized_GNEPc_gurobi_MA(p);
%         [p,o]= solveCentralized_GNEPc_gurobi(p)
%        [p,o]= solveCentralized_GNEPc(p);
        
        %% Stage 2
        o = solve_binary(o,p);

        [o,e] = solve_pressure3(o,p);

        o.pen_w(r) = gamma(r);
        er(r) = o.gfv_max;
    %    n_iter(r) = o.n_iter;


        %% Update gamma
        [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);
        q.time(r+1) = toc;
        
        % compute cost
        for i=1:p.n
            [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
        end
        q.Jt(r) = sum(q.J(:,r));
        q.Pt(r) = sum(q.P(:,r));

        q.gamma = gamma;
        q.er = er;
        
        %q.nIter = n_iter;
        q.gammaLower = gammaLower;
        q.gammaUpper = gammaUpper;
        %save(['sim1_2802_res',num2str(rr)],'o','e','p_data','q')
        
        q.er_gf = gasFlow_error(p,o);
        %%
        if er(r) < 1e-5
            break
        end
    end
    %clearvars('o','e','q','p','p_data')
end
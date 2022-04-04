function p = initialize_GNEPc_MA(p)

%% Identify dimensions of decision variables (per time step, h)
    for i = 1:p.n
        p.nxc(i) = length(p.en.Nex{i});   % auxiliary consensus variables
        p.nx(i) = 5+p.en.noN(i)+p.nxc(i);
        p.nt(i) = p.gn.noN(i);
        p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
        p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
        p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);

    % 
        for jj = 1:p.gn.noN(i)
            p.id_phi{i}(jj) = p.nx(i) + 2+(2+p.r)*(jj-1)+1;
        end
    end
    p.id_dg = 1;
    p.id_ch = 2;
    p.id_dh = 3;
    p.id_gu = 4;
    p.id_th = 5;
    % y vars. must +p.nx(i)
    p.id_psi = 1;
    p.id_gs = 2;


    
       
    % Generate matrices for linear constraints       
    p = build_mat_iegs_opt(p); 
    
end
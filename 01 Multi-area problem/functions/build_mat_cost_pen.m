function p = build_mat_cost_pen(p)

for i = 1:p.n
%% Cost function 
    q = zeros(1,p.nx(i));
    q(p.id_dg) = p.en.q_dg(i);
    q(p.id_ch) = p.en.q_st(i);
    q(p.id_dh) = p.en.q_st(i);
    
    Q = diag(q);
    Qp = diag(q);
    
    Qh{i} = p.m.Sx{i}'*kron(eye(p.h),Q)*p.m.Sx{i};
    Qph{i} = p.m.Sx{i}'*kron(eye(p.h),Qp)*p.m.Sx{i};
    
    p.m.Qh{i} = sparse(Qh{i});
    p.m.Qph{i} = sparse(Qph{i});
    
    p.m.H{i} = p.Alpha{i} + 2*Qh{i};
    

    
    cc = 1;

    c = zeros(p.nu(i),1);
    c(p.id_dg) = p.en.c_dg(i);
    c(p.id_ch) = p.en.c_st(i);
    c(p.id_dh) = p.en.c_st(i);
    c(p.nx(i)+p.id_gs) = p.gn.c_gs(i);
    
    % only for alg of he'18, penalty-based approach. 
    if p.Gamma_pen_flag == 1
        eta = p.gamma_pen; % weight
        for j = 1:p.gn.noN(i)
            c(p.nx(i)+2+2*j) = eta;
        end
    end
    
    p.m.ch{i} = kron(ones(p.h,1),c) ;
    
    % infinity norm penalty
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        p.m.ch{i} = p.m.ch{i} + p.m.Spen{i,j}'*ones(p.h,1)*p.pen;
    end
end
end
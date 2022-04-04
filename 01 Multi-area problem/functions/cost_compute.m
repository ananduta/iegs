 function [Ji,Pi] = cost_compute(o,p,i)

%     % compute aggregative terms
%     o.sigma_mg = p.en.sumPd(1:p.h);
%     o.sigma_gu = p.gn.sumGd(1:p.h);
%     for j=1:p.n
%         o.sigma_mg = o.sigma_mg + o.p_mg{j};
%         o.sigma_gu = o.sigma_gu + o.d_gu{j};
%     end


    
    c1 = [];
%     cc = 1;
%     c = [p.en.c_dg(i) p.en.c_st(i) p.en.c_st(i) p.en.d_mg_l p.gn.d_gu_l 0 0 0 zeros(1,p.en.noN(i))]';
%     for jj=1:p.en.noN(i)
%         j = p.en.N{i}(jj);
%         c(4+cc,1) = p.en.c_tr(i,j);
%         cc = cc+1;
%     end
    c = zeros(p.nu(i),1);
    c(p.id_dg) = p.en.c_dg(i);
    c(p.id_ch) = p.en.c_st(i);
    c(p.id_dh) = p.en.c_st(i);
    c(p.nx(i)+p.id_gs) = p.gn.c_gs(i);

    
    
    
    Ji = o.u{i}'*p.m.Qh{i}*o.u{i} + kron(ones(p.h,1),c)'*o.u{i} +p.h*p.en.c2_dg(i);
    
        
    Pi = Ji;
end
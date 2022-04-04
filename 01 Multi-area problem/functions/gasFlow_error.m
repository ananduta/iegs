function er = gasFlow_error(p,o)
% Compute error of gas flow non-linear equations
% W. Ananduta
% 01/03/2022
c = 1;
for i = 1:p.n
    for jj = 1:length(p.gn.Nin{i})
        j = p.gn.Nin{i}(jj);
        
        for h = 1:p.h
        % gas flow equation
        
        if o.psi{i} >= o.psi{j}
            sig = -1;
        else
            sig = 1;
        end
        
        
        flow(c,h) = sqrt(abs(o.psi{i}(h)-o.psi{j}(h)));
        gf(c,h) = ((1/sqrt(p.gn.cf(i,j))*(o.phi{i,j}(h))) - sig*sqrt(abs(o.psi{i}(h)-o.psi{j}(h))))/flow(c,h);
        gf_abs(c,h) = (1/sqrt(p.gn.cf(i,j))*(o.phi{i,j}(h)) - sig*sqrt(abs(o.psi{i}(h)-o.psi{j}(h))));
        g_f(i,j,h) = gf(c,h);
        
        if abs(gf(c,h)) < 1e-5
            gf_tight(c,h) = 1;
        else
            gf_tight(c,h) = 0;
        end
        
        c = c+1;
        end
    end
    
    
end
er.max = max(max(abs(gf)));
er.mean = mean(mean(abs(gf)));
er.gf = gf;
er.gf_abs = gf_abs;
er.gf_tight = gf_tight;
er.flow = flow;
end



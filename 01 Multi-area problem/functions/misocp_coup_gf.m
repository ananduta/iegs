function misocpCoup = misocp_coup_gf(p)
% build matrices for linear coupling constraints of mixed-integer SOCP reformulation
% of the gas-flow equations
% W. Ananduta
% 16/02/2022
for i = 1:p.n
	for jj = 1:p.gn.noN(i)
        j=p.gn.N{i}(jj);
        
        id_psi = p.nx(i) + 1;
        id_gamma =  p.nx(i) + 2+(2)*jj;
        id_alpha =  p.nx(i) + p.ny(i) + jj;
        
        
        id_psi_j = p.nx(j) + 1;
        
        %(33)
        E1i = zeros(1,p.nu(i));
        E1i(1,id_gamma) = -1;
        E1i(1,id_psi) = 1;
        E1i(1,id_alpha) = 2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        
        E1j = zeros(1,p.nu(j));
        E1j(1,id_psi_j) = -1;
        
        f1 = 0;
    
        
        %(34)
        E2i = zeros(1,p.nu(i));
        E2i(1,id_gamma) = -1;
        E2i(1,id_psi) = -1;
        E2i(1,id_alpha) = 2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        
        E2j = zeros(1,p.nu(j));
        E2j(1,id_psi_j) = 1;
        
        f2 = 2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        %(35)
        E3i = zeros(1,p.nu(i));
        E3i(1,id_gamma) = 1;
        E3i(1,id_psi) = -1;
        E3i(1,id_alpha) = -2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        
        E3j = zeros(1,p.nu(j));
        E3j(1,id_psi_j) = 1;
        
        f3 = 0;
        
        %(36)
        E4i = zeros(1,p.nu(i));
        E4i(1,id_gamma) = 1;
        E4i(1,id_psi) = 1;
        E4i(1,id_alpha) = -2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        
        E4j = zeros(1,p.nu(j));
        E4j(1,id_psi_j) = -1;
        
        f4 = -2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        Eci = [E1i;E2i;E3i;E4i];
        Ecj = [E1j;E2j;E3j;E4j];
        fc = [f1;f2;f3;f4];
        
        Ec{i,j}{1} = kron(eye(p.h),Eci);
        Ec{i,j}{2} = kron(eye(p.h),Ecj);
        fc{i,j} = kron(ones(p,h,1),fc);
        
        
        
        
    end
end
misocpCoup.Ec = Ec;
misocpCoup.fc = fc;

end
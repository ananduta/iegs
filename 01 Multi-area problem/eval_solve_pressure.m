% relationship between particular sol. of solve_pressure and flow solution.
% 26/10/21

clear all
clc
%close all

n = 30; gn.n = n;
sp = 0.7;
A = randconG(n,sp);
G = graph(A);
G = minspantree(G);
gn.Adj = adjacency(G);
[gn.N,gn.noN] = id_neigh(gn.Adj);
normA=norm(full(gn.Adj));
maxlamA = max(eig(full(gn.Adj)));
gn.r = 10;
gn.phi_max = (30+(20*rand))*gn.n*ones(gn.n);
gn.cf = (35+(15*rand))*ones(n);

sam = 0.01;
par = sam:sam:1;


%%
c = 1;
max_c = 10;
for k = 1:length(par)
    for l = 1:max_c
        
        phi{c} = zeros(n);
        
        for i = 1:n
            for jj = 1:gn.noN(i)
            j = gn.N{i}(jj);
                if j>i
                    phi_min = (par(k)-sam)*gn.phi_max(i,j);
                    phi_max = par(k)*gn.phi_max(i,j);
                    phi{c}(i,j) = phi_min + (phi_max-phi_min)*rand;
                end
            end
        end
        
        
        phi{c} = phi{c} - phi{c}';
        c = c+1;
    end
end

e = 1;
for c = 1:length(phi)
    for i=1:n
        for jj = 1:gn.noN(i)
            j = gn.N{i}(jj);


            % define the parameters of the affine function at each region
            nf = @(y) y^2/gn.cf(i,j);
            pwaf = pwa_approx_nf(gn.r,-gn.phi_max(i,j),gn.phi_max(i,j),nf);
            
            % \delta_(i,j)^\psi
            if phi{c}(i,j) <= 0
                delta_psi_s(i,j) = 1;
            else
                delta_psi_s(i,j) = 0;
            end


            % \delta_(i,j)^m            

            for m = 1:gn.r

                if phi{c}(i,j) >= pwaf.m(m) && phi{c}(i,j) <= pwaf.M(m)
                    delta_s(i,j,m) = 1;
                else
                    delta_s(i,j,m) = 0;
                end
            end
            
            b_ij = 0;
            for m = 1:gn.r
                b_ij = b_ij - delta_s(i,j,m)*(pwaf.a(m)*phi{c}(i,j) + pwaf.b(m));
            end

            A_ij = diag(2*delta_psi_s(i,j)-1);
            
            E(e,i) = A_ij;
            E(e,j) = -A_ij;
            b(e,1) = b_ij;
            
            e = e+1;
        end
    end
    
    pre2(c,1) = norm(phi{c}(:),2);
    x02(c,1) = norm(pinv(E)*b,2);
    
    pre_i(c,1) = norm(phi{c}(:),inf);
    x0_i(c,1) = norm(pinv(E)*b,inf);
end

figure
scatter(pre2,x02)
xlabel('$|| \phi||_2$','Interpreter','latex')
ylabel('$|| x_0||_2$','Interpreter','latex')
figure
scatter(pre_i,x0_i)
xlabel('$|| \phi||_\infty$','Interpreter','latex')
ylabel('$|| x_0||_\infty$','Interpreter','latex')
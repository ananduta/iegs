% IEGS Network generation
% random gas network

% W. Ananduta
% 18/06/2021

% Generate an integrated electrical and gas (IEGS) network
% Set the parameters of the network
% Set the gas and electrical loads

% time horizon
p.h = 1;
t_start = 12;

%
n_agents = 33;
p.n = n_agents;

% scaling
p.scale = 1;

% load network data
load('case33b20n.mat')

%% Electrical network (en)
% set time horizon
np.h= p.h;

np.n = p.n;



% Set B and G

np.Bnet = zeros(np.n);
np.Gnet = zeros(np.n);
np.Adj = zeros(np.n);


tab = table2array(linesE);

for i= 1:length(tab(:,1))
    np.Adj(tab(i,2),tab(i,3)) = 1;
    
    z = tab(i,4) + tab(i,5)*sqrt(-1);
    
    
    np.Bnet(tab(i,2),tab(i,3)) = abs(imag(1/z));
    np.Gnet(tab(i,2),tab(i,3)) = abs(real(1/z));
end
np.Adj = np.Adj' + np.Adj;
[np.N,np.noN] = id_neigh(np.Adj);

np.v_op = 12.66; 

np.scaling = 1e6;      %    EDIT HERE
p.scaling = np.scaling;
np.Bnet = sparse(np.Bnet + np.Bnet');
np.Bnet = np.v_op^2*1000*np.Bnet/np.scaling;
np.Gnet = sparse(np.Gnet + np.Gnet');
np.Gnet = np.v_op*np.Gnet/np.scaling;


np.pas_ag = 0;

% Define partition
p.nA = 3;
p.Ae{1} = [[1:1:5],[19:1:25]];
p.Ae{2} = [6:1:18];
p.Ae{3} = [26:1:33];

p.Ag{1} = [1:1:7];
p.Ag{2} = [8:1:16];
p.Ag{3} = [17:1:20];

p.A_flag = 1;
% Interconnection between electrical and gas networks


gn.no_nodes = 20;
%sp = sp_set;

gn.scale = p.scale;

% Randomly pairing a bus in the electrical network and a node in the gas
% network (no nodes <= no busses)
pairList = zeros(1,p.n);
for i = 1:p.nA
    no_nodes = length(p.Ag{i});
    pList_id = [1:no_nodes zeros(1,length(p.Ae{i})-no_nodes)];
    pList_id = pList_id(randperm(length(pList_id)));
    
    for j = 1:length(p.Ae{i})
        if pList_id(j) > 0
            pairList(p.Ae{i}(j)) = p.Ag{i}(pList_id(j)); 
            pList_A{i}(j) = p.Ag{i}(pList_id(j));
        end
    end

end
gn.pairList = pairList;
p.pList_A = pList_A;
p.pairList = pairList;

% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_et = 1;
 ag_et = 1; 
 np.et = zeros(np.n,1);
 for i = 1:n_et
     np.et(ag_et(i)) = 1;
 end
 
 
% assign components to each agent
n = np.n;
b = np.pas_ag;
np.b = b;
% type of load
% randomly assign the type of load profile   
np.t_lpr = randi([3 6], n+b,1);
% 0 = no load
%np.t_lpr(n+1) = 0;
%np.t_lpr(n+7) = 0;
%np.t_lpr(n+8) = 0;
    
% storage units (no storage units in [Sousa, et. al.,2019] and [Le Cadre, et al, 2019]
%np.st_un = zeros(n,1);
n_st=floor(np.n/3);
 ag_st = randperm(np.n,n_st); 
 np.st_un = zeros(n,1);
 for i = 1:n_st
     np.st_un(ag_st(i)) = 1;
 end

% subset of gas-powered dispatchable units
n_gu = 5;
busnode = find(pairList~=0);

ag_gu_id = randperm(length(busnode),n_gu);
ag_gu = zeros(n_gu,1);
for i = 1:n_gu
    ag_gu(i) = busnode(ag_gu_id(i));
end
np.dgu_un = zeros(n,1);
for i = 1:length(ag_gu)
    np.dgu_un(ag_gu(i)) = 1;
end
 
% dispatchable units
%np.d_un = randi([0 1], n,1);
%np.d_un = [ones(1, 10) zeros(1,n-10)];

 
 %ag_dg = [2,4,7,8,14,15,18,24,25,32];
 n_ngu = 5;
 busNotNode = find(pairList==0);
 ag_ngu_id = randperm(length(busNotNode),n_ngu);
 ag_ngu = zeros(n_ngu,1);
 for i = 1:n_ngu
     ag_ngu(i) = busNotNode(ag_ngu_id(i));
 end
 ag_dg = [ag_gu;ag_ngu];
 
 np.d_un = zeros(n,1);
 for i = 1:length(ag_dg)
     np.d_un(ag_dg(i)) = 1;
end
%np.d_un = zeros(n,1);



% PV generation units
np.r_un = zeros(n+b,1);
for i=1:n
    if np.st_un == 1
        np.r_un(i) = randi([0 4]);
    end
end



% assign parameters in the local constraints
np.scale=p.scale;
np = gen_param(np,ty); 



% generate power load and non-dispatchable profiles



tab = table2array(loadsE);

for i=1:length(tab(:,1))
    Pl = gen_load(np.t_lpr(i));
    loadRatioE = Pl/max(Pl);
    np.Pd(1:p.h,i) = tab(i,2)*loadRatioE(t_start:t_start+p.h-1);
end
%np.Pd(1:p.h,1) = zeros(p.h,1);   % zero load on the pcc 

np.Pd = np.Pd';
np.sumPd = sum(np.Pd(np.n+1:end,:))';



% initial condition of variables
np.init = 0;

% assign per-unit costs
np = gen_cost(np,tc); 

p.en = np;

clearvars('np');




%% Gas network (gn)
gn.n = p.n; 

% Adjacency matrix of gas network
maxBranch = 3;
gn.Adj = zeros(p.n);
cf = zeros(p.n);


tab = table2array(linesG);
for i = 1:length(tab(:,1))
    node1 = find(pairList==tab(i,2));
    node2 = find(pairList==tab(i,3));
    
    gn.Adj(node1,node2) = 1;
    cf(node1,node2) = tab(i,4);
end
cf = cf/gn.scale*sqrt(gn.scale);
cf = cf + cf';

gn.Adj = gn.Adj' + gn.Adj;

[gn.N,gn.noN] = id_neigh(gn.Adj);



gn.cf = cf.^2; %squared c^f




gn.Ngu = ag_gu; % Set of nodes in Gas network with gas-powered units
gn.gfg = p.en.dgu_un; % Indicator whether the nodes have gas-powered units


% Set nodes connected to gas source
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_gs = 2;
node_gs = [1,8];
for i = 1:length(node_gs)
    ag_gs(i) = find(pairList==node_gs(i));
end
gn.gs = zeros(p.n,1);
 for i = 1:n_gs
     gn.gs(ag_gs(i)) = 1;
 end


% generate gas demand
gn.Gd_max = 6/gn.scale;

gn.Gdem = zeros(p.n,p.h);


tab = table2array(loadsG);
loadratioG = table2array(loadratioG);

for i = 1:length(tab(:,1))
    node = find(pairList==tab(i,2));
    gn.Gdem(node,:) = tab(i,3)/gn.scale*loadratioG(t_start:t_start+p.h-1,2);
   
end
%gn.Gdem(1:p.h,1) = zeros(p.h,1);
gn.Gdem = (0.55+0.15*rand)*gn.Gdem;
gn.sumGd = sum(gn.Gdem(gn.n+1:end,:))';

% Set constraint parameters (pc)
gn.pres = pressuresG;
gn = param_cons_gasNetwork_20n(gn);

% per-unit cost of gas
%gn.d_gu = 0.1412*0.05/mean(1+sum(gn.Gdem(gn.n+1:end,:)));
%gn.d_gu_l = 0.001;
gn.c_gs = 0.05*ones(gn.n,1);

% Set number of regions for PWA gas-flow model
p.r = 20;
gn.r = p.r;
p.gn = gn;


% Default parameters for some algorithms
p.pen = 0; % for our algorithm
p.Gamma_pen_flag = 0; % for he'18, if penalty-based used, then set to 1


clearvars('gn');





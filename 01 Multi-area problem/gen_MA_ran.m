%% Generate multi-area IEGS network
% W. Ananduta
% 21/02/2022


% Generate the network and its components, and define the parameters of the network and the components
% Also define partitions (areas) of the network
% run('gen_iegs_33b_20n_MA.m')
%run('gen_iegs_33b_20n.m')

%run('gen_iegs_12n_MA_ran.m')
%run('gen_iegs_73n_MA_ran.m')
run('gen_iegs_472n_MA_ran.m')
if p.A_flag == 0
    p.nA = 1;
    p.Ae{1} = [1:1:p.n];
    p.pList_A{1} = p.gn.pairList;
end

% index area of each bus and node
p.Ae_id = zeros(1,p.n);
p.Ag_id = zeros(1,p.gn.no_nodes);
for i = 1:p.n
    for j = 1:p.nA
        if ~isempty(find(p.Ae{j}==i))
            p.Ae_id(i) = j;
            
        end
        
        if ~isempty(find(p.pList_A{j}==i))
            p.Ag_id(i) = j;
        end
    end
end


% Identify internal and external neighbors of each bus and node
for i = 1:p.n
    c_in = 1;
    c_ex = 1;
    p.en.Nin{i} = [];
    p.en.Nex{i} = [];
    for jj = 1:p.en.noN(i)
        j = p.en.N{i}(jj);
        if p.Ae_id(j) == p.Ae_id(i)
            p.en.Nin{i}(c_in) = j;
            c_in = c_in + 1;
        else
            p.en.Nex{i}(c_ex) = j;
            c_ex = c_ex + 1;
        end
    end
    
    c_in = 1;
    c_ex = 1;
    p.gn.Nin{i} = [];
    p.gn.Nex{i} = [];
    i_id = p.gn.pairList(i); 
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        j_id = p.gn.pairList(j);
        if p.Ag_id(j_id) == p.Ag_id(i_id)
            p.gn.Nin{i}(c_in) = j;
            c_in = c_in + 1;
        else
            p.gn.Nex{i}(c_ex) = j;
            c_ex = c_ex + 1;
        end
    end
end 


 
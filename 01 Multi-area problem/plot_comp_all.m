%% plotting preliminary comparison results
clear all
% load data
ndata = 10;
for i = 1:ndata
    dat{i} = ['sim3_12n_0203_',num2str(i),'.mat'];
    %dat{i} = ['sim2_73n_0703_',num2str(i),'.mat'];
    %dat{i} = ['sim2_472n_0903_',num2str(i),'.mat'];
end

% collect data into useful arrays
Jt = zeros(ndata,3);
time = zeros(ndata,3);
gf = [];
id = [1,3,2];
for i = 1:ndata
    load(dat{i});
    
    % cost
    J(1,1) = q1{id(1)}.Jt(end);
    J(1,2) = q1{id(2)}.Jt(end);
    J(1,3) = q2{id(3)}.Jt(end);
    J = J/J(1,1);
    Jt(i,:) = J;
    
    % gas-flow violation
    gf = [gf;
          abs(q1{id(1)}.er_gf.mean) abs(q1{id(2)}.er_gf.mean) abs(q2{id(3)}.er_gf.mean)];
    
    % computational time
    time(i,1) = sum(q1{id(1)}.time);
    time(i,2) = sum(q1{id(2)}.time);
    time(i,3) = sum(q2{id(3)}.time);
end
%%

figure
subplot(3,1,1)
boxplot(Jt)
title('Cost (proportional)')
name = {'MISOCP';'MISOCP pen';'Our Method'};
set(gca,'xticklabel',name)

subplot(3,1,2)
boxplot(gf)
title('Average gas flow deviation')
name = {'MISOCP';'MISOCP pen';'Our Method'};
set(gca,'xticklabel',name)

subplot(3,1,3)
boxplot(time)
title('Computational time [s]')
name = {'MISOCP';'MISOCP pen';'Our Method'};
set(gca,'xticklabel',name)
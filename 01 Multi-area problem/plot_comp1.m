%% plotting preliminary comparison results
%clear all
% default data
%dat = 'sim_12n_0303b.mat';
%dat = 'sim_73n_0303a.mat';
%load(dat)
r = 1;
% Cost
cost(1,1) = q1{1}.Jt(end);
%cost(2,1) = q1{3}.Jt(end);
dat_pen = 1;
for i=1:r
    cost(i+dat_pen,1) = q{i}.Jt(end);
end
cost = cost/max(cost);

% max Gas flow violation
fv_max(1,1) = q1{1}.er_gf.max;
%fv_max(2,1) = q1{3}.er_gf.max;
for i=1:r
    fv_max(i+dat_pen,1) = q{i}.er_gf.max;
end

% mean Gas flow violation
fv_mean(1,1) = q1{1}.er_gf.mean;
%fv_mean(2,1) = q1{3}.er_gf.mean;
for i=1:r
    fv_mean(i+dat_pen,1) = q{i}.er_gf.mean;
end
% computational time
time(1,1) = sum(q1{1}.time);
%time(2,1) = sum(q1{3}.time);
for i=1:r
    time(i+dat_pen,1) = sum(q{i}.time);
end

figure
subplot(3,1,1)
bar(cost)
title('Cost (proportional)')
name = {'MISOCP pen';'r=100'};
set(gca,'xticklabel',name)

subplot(3,1,2)
bar(fv_mean)
title('Average flow violation (proportional)')
name = {'MISOCP';'r=100'};
set(gca,'xticklabel',name)

subplot(3,1,3)
bar(time)
title('Computational time [s]')
name = {'MISOCP';'r=100'};
set(gca,'xticklabel',name)
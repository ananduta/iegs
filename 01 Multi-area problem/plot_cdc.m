%% plotting preliminary comparison results
clear all
%close all
% load data
ndata = 100;
%load('sim5_73n_1803.mat')
%load('sim3_73n_1103.mat')
load('sim5_472n_1903.mat')
% collect data into useful arrays
Jt = zeros(ndata,3);
time = zeros(ndata,3);
gf = [];
id = 2;
for i = 1:ndata
  
    
    % cost
    J(1,1) = q1{i}.Jt(end);
    J(1,2) = q2{i}.Jt(end);
    J(1,3) = q3{i,id}.Jt(end);
%     J(1,3) = q3{i}.Jt(end);
    if J(1,1) > 1e10
        J = ones(1,3);
    else
        J = J/J(1,1);
    end
    
    Jt(i,:) = J;
    
    % gas-flow violation
    gf = [gf;
         abs(q1{i}.er_gf.mean) abs(q2{i}.er_gf.mean) abs(q3{i,id}.er_gf.mean)];
%     gf = [gf;
%           abs(q1{i}.er_gf.mean) abs(q2{i}.er_gf.mean) abs(q3{i}.er_gf.mean)];
    % computational time
    time(i,1) = sum(q1{i}.time);
    time(i,2) = sum(q2{i}.time);
    time(i,3) = sum(q3{i,id}.time);
%     time(i,3) = sum(q3{i}.time);
end
set(groot,'defaultAxesFontSize',14)
% set(groot,'defaultAxesFont','Times')
 set(groot,'defaultAxesTickLabelInterpreter','latex');
%set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

figure
subplot(3,1,1)
hold on, grid on, box on,

h1 = bar(mean(Jt))
set(h1,'FaceColor',[0, 0.75, 0.75])
h=boxplot(Jt);
set(h,{'linew'},{1.5})
title('Cost (proportional)','Interpreter','Latex','FontSize',14)
name = {'MISOC';'pen-MISOC';'Alg. 1'};
set(gca,'xticklabel',name)

subplot(3,1,2)
hold on, grid on, box on,

h1 = bar(mean(gf))
set(h1,'FaceColor',[0, 0.75, 0.75])
h = boxplot(gf);
set(h,{'linew'},{1.5})
title('Average gas flow deviation','Interpreter','Latex','FontSize',14)
name = {'MISOC';'pen-MISOC';'Alg. 1'};
set(gca,'xticklabel',name)

subplot(3,1,3)
hold on, grid on, box on,

h1 = bar(mean(time))
set(h1,'FaceColor',[0, 0.75, 0.75])
h=boxplot(time);
set(h,{'linew'},{1.5})
title('Computational time [s]','Interpreter','Latex','FontSize',14)
name = {'MISOC';'pen-MISOC';'Alg. 1'};
set(gca,'xticklabel',name)

%% second plot evaluation of error vs #pwa regions
%load('sim5_472n_1903.mat')
load('sim5_73n_1803.mat')
gf = [];
tg = [];
for i = 1:ndata
    
    gf = [gf;
         abs(q3{i,1}.er_gf.mean) abs(q3{i,2}.er_gf.mean) abs(q3{i,3}.er_gf.mean) abs(q3{i,4}.er_gf.mean)];
    
    
end

figure
%subplot(2,1,1)
hold on, grid on, box on,
h1 = bar(mean(gf))
set(h1,'FaceColor',[0, 0.75, 0.75])
h = boxplot(gf);
set(h,{'linew'},{1.5})
title('Average gas flow deviation','Interpreter','Latex','FontSize',14)
name = {'20';'40';'60';'80'};
set(gca,'xticklabel',name)
xlabel('Number of PWA regions','Interpreter','Latex','FontSize',14)
% load('sim5_472n_1903.mat')
% gf = [];
% tg = [];
% for i = 1:ndata
%     gf = [gf;
%          abs(q3{i,1}.er_gf.mean) abs(q3{i,2}.er_gf.mean) abs(q3{i,3}.er_gf.mean) abs(q3{i,4}.er_gf.mean)];
%     
%     
% end
% 
% subplot(2,1,2)
% hold on, grid on, box on,
% h1 = bar(mean(gf))
% set(h1,'FaceColor',[0, 0.75, 0.75])
% h = boxplot(gf);
% set(h,{'linew'},{1.5})
% title('472-node-4-area test case','Interpreter','Latex','FontSize',14)
% name = {'20';'40';'60';'80'};
% set(gca,'xticklabel',name)




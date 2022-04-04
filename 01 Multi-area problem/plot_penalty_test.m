%% plot analysis of penalty-based method
clear all
clc
close all
data{1} = 'sim0911.mat';

pen = [0,1,10,50,100,120,125,130,135,140,160,180,200,300,400,500,1000];
for  d = 2:length(pen)
    data{d} = ['sim_0911_',num2str(pen(d)),'.mat'];
end


for d = 1:length(data)
    load(data{d})
    phi_n(d) = (norm(o.phi{5,6},inf));
    err(d) = o.gfv_max;
    
    for i=1:p.n
        J(i,d) = cost_compute(o,p,i);
    end
    
end

figure
%grid on, box on
subplot(2,1,1)
grid on
box on
plot(pen,phi_n,'o-','LineWidth',2)
xlabel('penalty coef','Interpreter','Latex','FontSize',14)
ylabel('$||$flows$||_\infty$','Interpreter','Latex','FontSize',14)

subplot(2,1,2)
grid on
box on
plot(pen,err,'o-','LineWidth',2)
xlabel('penalty coef','Interpreter','Latex','FontSize',14)
ylabel('$||$error gas flow eq.$||_\infty$','Interpreter','Latex','FontSize',14)


figure
grid on, box on, hold on
for i =1:p.n
    plot(pen,J(i,:),'o-','LineWidth',1.3)
end
plot(pen,sum(J),'o-','LineWidth',2)
xlabel('penalty coef','Interpreter','Latex','FontSize',14)
ylabel('cost','Interpreter','Latex','FontSize',14)
set(gca, 'YScale', 'log')
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Total')
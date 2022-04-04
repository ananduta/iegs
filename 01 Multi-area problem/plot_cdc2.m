%%
clear all
load('sim5_73n_1803.mat')
ndata = 100;
id = 4;
for id=1:4
    for i=1:ndata
        er(i,id)= q3{i,id}.er;
        if abs(er(i,id)) <= 1e-6
            flag(i,id) = 1;
        else
            flag(i,id) = 0;
        end
    end
end
%max(er)
%sum(flag)

gf = [];
tg = [];
for i = 1:ndata
    if sum(flag(i,:))==4
        gf = [gf;
             abs(q3{i,1}.er_gf.mean) abs(q3{i,2}.er_gf.mean) abs(q3{i,3}.er_gf.mean) abs(q3{i,4}.er_gf.mean)];
    end
    time(i,1) = sum(q3{i,1}.time);
    time(i,2) = sum(q3{i,2}.time);
    time(i,3) = sum(q3{i,3}.time);
    time(i,4) = sum(q3{i,4}.time);
    
end

load('sim5_472n_1903.mat')
ndata = 100;
id = 4;
for id=1:4
    for i=1:ndata
        er(i,id)= q3{i,id}.er;
        if abs(er(i,id)) <= 1e-6
            flag1(i,id) = 1;
        else
            flag1(i,id) = 0;
        end
    end
end
%max(er)
%sum(flag)

gf1 = [];

for i = 1:ndata
    if sum(flag(i,:))==4
        gf1 = [gf1;
             abs(q3{i,1}.er_gf.mean) abs(q3{i,2}.er_gf.mean) abs(q3{i,3}.er_gf.mean) abs(q3{i,4}.er_gf.mean)];
    end
    time1(i,1) = sum(q3{i,1}.time);
    time1(i,2) = sum(q3{i,2}.time);
    time1(i,3) = sum(q3{i,3}.time);
    time1(i,4) = sum(q3{i,4}.time);
end

gf = [gf;gf1];
flag = [flag;flag1];
%sflag = sum([flag;flag])/2;
figure
subplot(3,1,1)
hold on, grid on, box on,
h1 = bar([20:20:80],sum(flag)/2)
set(h1,'FaceColor',[0, 0.75, 0.75])
title('Frequency of optimal solutions found','Interpreter','Latex','FontSize',14)
%name = {'20';'40';'60';'80'};
%set(gca,'xticklabel',name)
xlabel('Number of PWA regions','Interpreter','Latex','FontSize',14)
ylabel('$\%$','Interpreter','Latex','FontSize',14) 
xlim([10,90])
xticks([20:20:80])
yticks([0:20:100])
subplot(3,1,2)
hold on, grid on, box on,
h1 = bar(mean(gf))
set(h1,'FaceColor',[0, 0.75, 0.75])
h = boxplot(gf);
set(h,{'linew'},{1.5})
title('Average gas flow deviation','Interpreter','Latex','FontSize',14)
name = {'20';'40';'60';'80'};
set(gca,'xticklabel',name)


subplot(3,1,3)
hold on, grid on, box on,

bar([20:20:80],[mean(time); mean(time1)])
%set(h1,'FaceColor',[0, 0.75, 0.75])
%h=boxplot(time);
%set(h,{'linew'},{1.5})
title('Computational time [s]','Interpreter','Latex','FontSize',14)
%name = {'20';'40';'60';'80'};
set(gca,'xticklabel',name)
xlabel('Number of PWA regions','Interpreter','Latex','FontSize',14)
legend('medium case','large case','Interpreter','Latex','FontSize',13)
xlim([10,90])
xticks([20:20:80])
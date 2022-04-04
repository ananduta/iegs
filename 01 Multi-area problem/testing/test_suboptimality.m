% test computation suboptimality
% W. Ananduta

% clear all
% clc

%load('sim_0911_1.mat')
%load('sim_0911_140.mat')

load('param_case.mat')
V = zeros(p.n,1);
Vc = zeros(p.n,1);
Ji_n = zeros(p.n,1);
for i = 1:p.n
    [V(i),Vc(i),Ji_n(i)] = NI_function(o,p,i);
end


function [xi,yij, sn, k, Obj]  = optimiza_CRSM(EmProba, DC, Sbar, His, b, tau)



%% 由于每个cluster中的需求点数不同，因此只能手动添加约束，先以K=3为例进行求解。
% Define the variables
lamda = sdpvar(1,1,'full');
xi = binvar(size(DC,1),1,'full');
yij1 = binvar(size(DC,1),size(DiffertKK{1,2},1),'full'); %
yij2 = binvar(size(DC,1),size(DiffertKK{2,2},1),'full'); %
yij3 = binvar(size(DC,1),size(DiffertKK{3,2},1),'full'); %
% yij4 = binvar(size(DC,1),size(DiffertKK{4,2},1),'full'); %
% yij5 = binvar(size(DC,1),size(DiffertKK{5,2},1),'full'); %
alpha1 = sdpvar(size(DiffertKK{1,1},1), 1,'full');
alpha2 = sdpvar(size(DiffertKK{2,1},1), 1,'full');
alpha3 = sdpvar(size(DiffertKK{3,1},1), 1,'full');
% alpha4 = sdpvar(size(DiffertKK{4,1},1), 1,'full');
% alpha5 = sdpvar(size(DiffertKK{5,1},1), 1,'full');

%% Objective
Obj = 0;
F = [];  % Initialize the costraints





end


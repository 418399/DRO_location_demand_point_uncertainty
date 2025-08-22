clc;clear;
tic

load His ;
%   load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
 [His] = Distance(His) ;
b = 10;   % The max number of opened DC
aa = 29;  % number of demand points in x-axes
bb = 29  ;% ...in y-axes
MultiFactor = 1.4;
aaa = 7 ;
bbb = 4 ;
K=3;

minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
% minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
% maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
lengthY = maxlo(2) - minlo(2) ;

Disc_X = [minlo(1):lengthX/aa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_Y = [minlo(2):lengthY/bb:maxlo(2)] ;

LC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_X,2)
    for j = 1:size(Disc_Y,2)
        LC = [LC ; Disc_X(i) , Disc_Y(j) ] ;
    end
end

% DC
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end
%[DC] = Distance(DC) ;
EmProba = 1/size(His,1) ;
% 得到 Z_0
[z_0 x_0] = getZ0(EmProba, DC, His, b) ;
% 选择 tau  +++++++++++++++++++++++++

tau = z_0 * MultiFactor ;
% No_initialLC = round(0.2* size(LC,1) ) ;
% random_Num = randperm(size(LC,1),No_initialLC) ;
%  [idx,Center] = kmeans(His,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster
load idx;
load Center ;
His = [His idx ] ;  % add the cluster number of each points in the 3th column

%% Judge the candidate points belongs to which sub-region
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
Vcells = Voronoi(Center ,area)  ;
hold on
plot(His(:,1),His(:,2),'v')
%
close()
temp = 0;
for k = 1:K
    [in,on] = inpolygon(LC(:,1), LC(:,2), Vcells{k,1},Vcells{k,2});
    in_poly = double(in)  ;
    temp = temp + in_poly*k ;  % convert to the NO. of each demand point
end
LC(:,3) = temp ; % add the NO. to the 3th column

%% Didstinguish the clusters in different sub-regions
for k = 1:K
    DiffertKK(k,1) = { His( find(His(:,3) == k), 1:2 )  }  ;  % His
    DiffertKK(k,2) = {  LC( find( LC(:,3) == k ), 1:2  ) }  ; % LC
end

%% 由于每个cluster中的需求点数不同，因此只能手动添加约束，先以K=3为例进行求解。
% Define the variables

xi = binvar(size(DC,1),1,'full');
yij1 = binvar(size(DC,1),size(DiffertKK{1,2},1),'full'); %
yij2 = binvar(size(DC,1),size(DiffertKK{2,2},1),'full'); %
yij3 = binvar(size(DC,1),size(DiffertKK{3,2},1),'full'); %
%   yij4 = binvar(size(DC,1),size(DiffertKK{4,2},1),'full'); %
%   yij5 = binvar(size(DC,1),size(DiffertKK{5,2},1),'full'); %
%
alpha1n = sdpvar(size(DiffertKK{1,1},1), 1,'full');
alpha2n = sdpvar(size(DiffertKK{2,1},1), 1,'full');
alpha3n = sdpvar(size(DiffertKK{3,1},1), 1,'full');
%   alpha4n = sdpvar(size(DiffertKK{4,1},1), 1,'full');
%   alpha5n = sdpvar(size(DiffertKK{5,1},1), 1,'full');

r = sdpvar(1,1,'full');

%% Objective
alpha_N_k_temp = 0 ; % 记录每个k循环下的第一个约束的temp
Obj = r;
F = [];  % Initialize the costraints
F = [F, sum(xi) == b] ;

for k = 1:K
    pk(k) = size(DiffertKK{k,1},1) / size(His,1) ;
    pro_nk(k) = 1/size(DiffertKK{k,1},1) ;  % the probability of N_k
    if k==1 % Constraints
        for n = 1:size(DiffertKK{k,1},1)
            alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha1n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij1(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij1(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha1n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij1(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    
    if k==2 % Constraints
        for n = 1:size(DiffertKK{k,1},1)
            alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha2n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij2(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij2(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha2n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij2(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    
    if k==3 % Constraints
        for n = 1:size(DiffertKK{k,1},1)
            alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha3n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij3(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij3(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha3n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij3(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==4 % Constraints
        for n = 1:size(DiffertKK{k,1},1)
            alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha4n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij4(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij4(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha4n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij4(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==5 % Constraints
        for n = 1:size(DiffertKK{k,1},1)
            alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha5n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij5(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij5(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha5n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij5(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    k
end

F = [F, alpha_N_k_temp <= tau];
F = [F, r >=0];
ops = sdpsettings('solver','cplex','verbose',3 );
%ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,Obj,ops );
solmp.info
x = value(xi);
r = value(r);


%% Verifing the location decisions via the generated test cases
load Test_points_55;
 %%Test results
Test_poi = [];
for ii = 1:10
    Test_poi = [Test_poi;  Test_points(:,2*ii-1:2*ii) ] ;
end

for i=1:size(DC,1)  % DC到所有DP的2范数
    for j = 1:size(Test_poi,1)
        save_2_norm(i,j) = norm( DC(i,:) -Test_poi(j,:) ,2)  ;
    end
end

temp_save_2_norm = save_2_norm( find(x==1),:)  ;
[min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % 对于每个z，找到与其最近的DC及距离
Save_results = min_value ;
Save_CVaR_Result = CVaRR(Test_poi, DC, x)  ;

% Print results
MAX_test = max(Save_results) ;
MEAN_test = mean(Save_results);
Save_CVaR_Result;
Print_results = [MAX_test ;  MEAN_test ; Save_CVaR_Result] 













% sx = [x;r];
% savex(ii) = {sx} ;
% end
% save('CRS_No_facility.mat','savex')
%alpha1n = value(alpha1n)
% hold on
%  %  gscatter(His(:,1),His(:,2),idx )
% plot(His(:,1),His(:,2),'v')
% hold on
% plot(DC(:,1),DC(:,2),'o')
toc

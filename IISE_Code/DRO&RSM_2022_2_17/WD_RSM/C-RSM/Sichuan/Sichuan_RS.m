clc;clear;tic
load Sichuan_data;
load Sichuan_DC;
%Sichuan_data = Distance(Sichuan_data);
b = 10;   % The max number of opened DC
aa = 39;  % number of demand points in x-axes
bb = 39  ;% ...in y-axes
K=  3 ;
% MultiFactor = 1.05;
MultiFactor = [1.2 1.3 1.4 1.5 1.6 ];

minlo = [min([Sichuan_data(:,1) ]) min([Sichuan_data(:,2) ])] ;
maxlo = [max([Sichuan_data(:,1) ]) max([Sichuan_data(:,2)])];
% minlo = [min([Sichuan_data(:,1) ; Sichuan_DC(:,1)]) min([Sichuan_data(:,2) ;Sichuan_DC(:,2)])] ;
% maxlo = [max([Sichuan_data(:,1) ;Sichuan_DC(:,1)]) max([Sichuan_data(:,2);Sichuan_DC(:,2)])];
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

EmProba = 1/size(Sichuan_data,1) ;  %每个历史点的发生概率
% 得到 Z_0
[z_0 x_0] = getZ0_distance(EmProba, Sichuan_DC, Sichuan_data, b) ;


No_initialLC = round(0.2* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;
%  [idx,Center] = kmeans(Sichuan_data,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster
load idxsichuan5;
load Centersichuan5 ;
Sichuan_data = [Sichuan_data idx ] ;  % add the cluster number of each points in the 3th column
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
Vcells = Voronoi(Center ,area)  ;
hold on
plot(Sichuan_data(:,1),Sichuan_data(:,2),'v')
% axis equal
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
    DiffertKK(k,1) = { Sichuan_data( find(Sichuan_data(:,3) == k), 1:2 )  }  ;  % His
    DiffertKK(k,2) = {  LC( find( LC(:,3) == k ), 1:2  ) }  ; % LC
end

for kkk = 1:1
    tau = z_0 * MultiFactor(kkk) ;
    %% 由于每个cluster中的需求点数不同，因此只能手动添加约束，先以K=3为例进行求解。
    % Define the variables
    xi = binvar(size(Sichuan_DC,1),1,'full');
    yij1 = binvar(size(Sichuan_DC,1),size(DiffertKK{1,2},1),'full'); %
%     yij2 = binvar(size(Sichuan_DC,1),size(DiffertKK{2,2},1),'full'); %
%     yij3 = binvar(size(Sichuan_DC,1),size(DiffertKK{3,2},1),'full'); %
%     yij4 = binvar(size(Sichuan_DC,1),size(DiffertKK{4,2},1),'full'); %
%     yij5 = binvar(size(Sichuan_DC,1),size(DiffertKK{5,2},1),'full'); %
%     yij6 = binvar(size(Sichuan_DC,1),size(DiffertKK{6,2},1),'full'); %
%     yij7 = binvar(size(Sichuan_DC,1),size(DiffertKK{7,2},1),'full'); %
    %yij8 = binvar(size(Sichuan_DC,1),size(DiffertKK{8,2},1),'full'); %
    
    %
    alpha1n = sdpvar(size(DiffertKK{1,1},1), 1,'full');
%     alpha2n = sdpvar(size(DiffertKK{2,1},1), 1,'full');
%     alpha3n = sdpvar(size(DiffertKK{3,1},1), 1,'full');
%     alpha4n = sdpvar(size(DiffertKK{4,1},1), 1,'full');
%     alpha5n = sdpvar(size(DiffertKK{5,1},1), 1,'full');
%     alpha6n = sdpvar(size(DiffertKK{6,1},1), 1,'full');
%     alpha7n = sdpvar(size(DiffertKK{7,1},1), 1,'full');
    %alpha8n = sdpvar(size(DiffertKK{8,1},1), 1,'full');
    r = sdpvar(1,1,'full');
    
    %% Objective
    alpha_N_k_temp = 0 ; % 记录每个k循环下的第一个约束的temp
    Obj = r;
    F = [];  % Initialize the costraints
    F = [F, sum(xi) == b] ;
    
    for k = 1:K
        pk(k) = size(DiffertKK{k,1},1) / size(Sichuan_data,1) ;
        pro_nk(k) = 1/size(DiffertKK{k,1},1) ;  % the probability of N_k
        if k==1 % Constraints
            for n = 1:size(DiffertKK{k,1},1)
                alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha1n(n) ;
                for j = 1:size(DiffertKK{k,2},1)
                    temp = 0;
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij1(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
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
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij2(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
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
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij3(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
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
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij4(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
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
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij5(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                        F = [F, yij5(i,j) <= xi(i)  ];  % allocation decision
                    end
                    F = [F, alpha5n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                    F = [F, sum(yij5(:,j)) == 1] ;  % allocation decision
                end
            end
        end

        if k==6 % Constraints
            for n = 1:size(DiffertKK{k,1},1)
                alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha6n(n) ;
                for j = 1:size(DiffertKK{k,2},1)
                    temp = 0;
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij6(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                        F = [F, yij6(i,j) <= xi(i)  ];  % allocation decision
                    end
                    F = [F, alpha6n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                    F = [F, sum(yij6(:,j)) == 1] ;  % allocation decision
                end
            end
        end   
        
        
        if k==7 % Constraints
            for n = 1:size(DiffertKK{k,1},1)
                alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha7n(n) ;
                for j = 1:size(DiffertKK{k,2},1)
                    temp = 0;
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij7(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                        F = [F, yij7(i,j) <= xi(i)  ];  % allocation decision
                    end
                    F = [F, alpha7n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                    F = [F, sum(yij7(:,j)) == 1] ;  % allocation decision
                end
            end
        end
        
        
        if k==8 % Constraints
            for n = 1:size(DiffertKK{k,1},1)
                alpha_N_k_temp = alpha_N_k_temp + pro_nk(k)*alpha8n(n) ;
                for j = 1:size(DiffertKK{k,2},1)
                    temp = 0;
                    for i = 1:size(Sichuan_DC,1)
                        temp = temp + yij8(i,j) * norm( Sichuan_DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                        F = [F, yij8(i,j) <= xi(i)  ];  % allocation decision
                    end
                    F = [F, alpha8n(n) >= pk(k) * ( temp - r*norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ) ]  ;   % every j  , k,  n has a corresponding constraint
                    F = [F, sum(yij8(:,j)) == 1] ;  % allocation decision
                end
            end
        end
        
        k
    end
    
    F = [F, alpha_N_k_temp <= tau];
    F = [F, r >=0];
    ops = sdpsettings('solver','cplex','verbose',2 );
    %ops.cplex.benders.strategy = 3 ;
    solmp=solvesdp(F,Obj,ops );
    solmp.info;
    
    x = value(xi)
    r = value(r)
    [averagedis] = sichuantest(x,Sichuan_DC)
    
    abc = toc
    save_res(:,kkk) = [x; r; averagedis; abc];
    
    
end




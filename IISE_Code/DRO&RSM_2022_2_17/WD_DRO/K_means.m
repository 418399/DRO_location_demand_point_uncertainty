clc;clear;
load His ;
load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
[DC] = Distance(DC) ;
b = 10;  % The max number of opened DC
aa = 29 ;  % number of demand points in x-axes
bb = 29 ;  % ...in y-axes
K = 3; % num. of sub-region
multiplier = 0.1 ;  %The factor to control the parameter \theta


[His] = Distance(His) ;

%% Generate candidate demand points
minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
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

%% DC
aaa = 7 ;
bbb = 4;
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end

%  get the max theta allover the region
EmProba = 1/size(His,1) ;
a = [];  % Find the MAX theta
for j = 1:size(LC,1)
    for n = 1:size(His,1)
        a(j,n) = norm(His(n,:) - LC(j,:),1  ) ;
    end
end
a = sum(a,2) ;
maxtheta = max(a) * EmProba;
%++++++++++++++++++++++++++++++++++++++++++++++++%

theta =multiplier * maxtheta ;
%+++++++++++++++++++++++++++++++++++++++++++++++%
   [idx,Center] = kmeans(His,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster
% load idx;
% load Center ;
His = [His idx] ;  % add the cluster number of each points in the 3th column

%% Judge the candidate points belongs to which sub-region
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
Vcells = Voronoi(Center,area)  ;
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
%% Model

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
for k = 1:K
    if k==1
        for n = 1:size(DiffertKK{k,1},1)
            Obj = Obj + EmProba * alpha1(n) ;
            % constraints: k==1
            for j = 1:size(DiffertKK{k,2},1) % Constraints  when k==1  yij(i,j) * norm(DC(i,1:2) - LC(j,1:2) ,2) ;
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij1(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij1(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha1(n) >= temp - lamda * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1)  ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij1(:,j)) == 1] ;  % allocation decision
            end  % j = J
        end   %n
    end  %k
    %
    if k==2
        for n = 1:size(DiffertKK{k,1},1)
            Obj = Obj + EmProba * alpha2(n) ;
            % constraints: k==2
            for j = 1:size(DiffertKK{k,2},1) % Constraints  when k==1  yij(i,j) * norm(DC(i,1:2) - LC(j,1:2) ,2) ;
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij2(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2)  ;% calculate the distance between DC and demand points
                    F = [F, yij2(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha2(n) >= temp - lamda * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1)  ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij2(:,j)) == 1] ;
            end
        end
    end
    %
    if k==3
        for n = 1:size(DiffertKK{k,1},1)
            Obj = Obj + EmProba * alpha3(n) ;
            % constraints: k==3
            for j = 1:size(DiffertKK{k,2},1) % Constraints  when k==1
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij3(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2)  ;% calculate the distance between DC and demand points
                    F = [F, yij3(i,j) <= xi(i)  ];
                end
                F = [F, alpha3(n) >= temp - lamda * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1)  ]  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij3(:,j)) == 1] ;
            end
        end
    end
    
    k
end  %  The for---end of  k = 1:K

%Obj = Obj * EmProba ;
Obj = Obj + theta * lamda;

% number of opened facilities
F = [F, sum(xi) == b] ;
F = [F, lamda >=0];
ops = sdpsettings('solver','gurobi','verbose',3 );
%ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,Obj,ops );
solmp.info;
Obj = value(Obj)
lamda = value(lamda);
yij1 =  value(yij1);
yij2 =  value(yij2);
yij3 =  value(yij3);
xi =  value(xi) ;

%% Plot the cluster results
% hold on
% plot(   LC(find(LC(:,3)==1),1),  LC(find(LC(:,3)==1),2),'p'   ) % plot the clustered descrite demand points
hold on
gscatter(His(:,1),His(:,2),idx )
hold on
plot(Center(:,1),Center(:,2),'kx')
%legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster Centroid')
hold on
plot(LC(:,1), LC(:,2), 'ko')
hold on
plot (DC(:,1), DC(:,2), 's') ;
hold on
open = find(xi ==1)  ;
plot( DC(open,1),DC(open,2) , 'ks','MarkerFaceColor','k') ;
axis equal
% 
% for i = 1:20
%     for j = 1:size(DiffertKK{1,2},1)
%         if yij1(i,j) == 1
%             plot([DC(i,1), DiffertKK{1,2}(j,1)], [DC(i,2), DiffertKK{1,2}(j,2)] ) ;
%             hold on
%         end
%     end
% end
% 
% 
% for i = 1:20
%     for j = 1:size(DiffertKK{2,2},1)
%         if yij2(i,j) == 1
%             plot([DC(i,1), DiffertKK{2,2}(j,1)], [DC(i,2), DiffertKK{2,2}(j,2)] ) ;
%             hold on
%         end
%     end
% end
% 
% 
% for i = 1:20
%     for j = 1:size(DiffertKK{3,2},1)
%         if yij3(i,j) == 1
%             plot([DC(i,1), DiffertKK{3,2}(j,1)], [DC(i,2), DiffertKK{3,2}(j,2)] ) ;
%             hold on
%         end
%     end
% end


%voronoi(Center(:,1), Center(:,2))
%axis equal
% set(gca,'XLim',[-97 -94]);%X轴的数据显示范围
% set(gca,'YLim',[27 31]);%Y轴的数据显示范围
%   title(['theta =' num2str(multi)])

   [Obj_Value1] = Func_Test_RSM0001(xi')
     [Obj_Value5] = Func_Test_RSM0005(xi')
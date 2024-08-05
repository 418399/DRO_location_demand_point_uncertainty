clc;clear;
tic
 
 
load His ;
%   load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
[His] = Distance(His) ;

b = 10;   % The max number of opened DC
aa = 39;  % number of demand points in x-axes
bb = 39  ;% ...in y-axes
MultiFactor = 1.1;
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
No_initialLC = round(0.2* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;
 %  [idx5,Center5] = kmeans(His,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster
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

%% 最重要的一个问题，由于每个subregion中的DC和demand points数量不一致，定义变量和约束时容易出错












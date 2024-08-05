 clc;clear;
tic
 
 
load Sichuan_data;
% %   load DC_20 ;
% Sichuan_data = Sourcedata1(1:20,:);  % Sichuan_datatorical data point
%[Sichuan_data] = Distance(Sichuan_data) ;

b = 10;   % The max number of opened DC
aa = 39;  % number of demand points in x-axes
bb = 39  ;% ...in y-axes
MultiFactor = 1.1;
  aaa = 7 ;
  bbb = 4 ;
K= 3 ;

minlo = [min([Sichuan_data(:,1) ]) min([Sichuan_data(:,2) ])] ;
maxlo = [max([Sichuan_data(:,1) ]) max([Sichuan_data(:,2)])];
% minlo = [min([Sichuan_data(:,1) ; DC(:,1)]) min([Sichuan_data(:,2) ;DC(:,2)])] ;
% maxlo = [max([Sichuan_data(:,1) ;DC(:,1)]) max([Sichuan_data(:,2);DC(:,2)])];
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
EmProba = 1/size(Sichuan_data,1) ;
% 得到 Z_0
[z_0 x_0] = getZ0(EmProba, DC, Sichuan_data, b) ;
% 选择 tau  +++++++++++++++++++++++++

tau = z_0 * MultiFactor ;
No_initialLC = round(0.2* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;
   [idx,Center] = kmeans(Sichuan_data,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster
%         load idx;
%         load Center ;
Sichuan_data = [Sichuan_data idx ] ;  % add the cluster number of each points in the 3th column

%% Judge the candidate points belongs to which sub-region
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
Vcells = Voronoi(Center ,area)  ;
hold on 
  plot(Sichuan_data(:,1),Sichuan_data(:,2),'v')
  
  
  
  
clc;clear;
tic
load His ;
load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
[DC] = Distance(DC) ;
b = 10;  % The max number of opened DC
aa = 19 ;  % number of demand points in x-axes
bb = 19 ;  % ...in y-axes
[His] = Distance(His) ;
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
bbb = 4 ;
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end


EmProba = 1/size(His,1) ;
%% get the max theta allover the region
a = [];  % Find the MAX theta
for j = 1:size(LC,1)
    for n = 1:size(His,1)
        a(j,n) = norm(His(n,:) - LC(j,:),1  ) ;
    end
end
a = sum(a,2) ;
maxtheta = max(a) * EmProba ; %  The max Wasserstein distance in the region
%+++++++++++++++++++++++++++++++++++++++++++
theta = 0.4* maxtheta ;
%+++++++++++++++++++++++++++++++++++++++++++
lamda = sdpvar(1,1,'full');
yij = binvar(size(DC,1),size(LC,1),'full'); %
xi = binvar(size(DC,1),1,'full');
sn = sdpvar(size(His,1),1,'full');

%% Objective
Obj = lamda * theta + EmProba*sum(sn) ;
%% Constraints
F = [];
% number of opened facilities
F = [F, sum(xi) == b] ;
% allocation decision
for j = 1:size(LC,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end

% Dual constraints
sumtemp = 0;
for j = 1:size(LC,1)
    for i = 1:size(DC,1)
        sumtemp = sumtemp + yij(i,j) * norm(DC(i,:) - LC(j,:),2) ;
    end
    for n = 1:size(His,1)
        F = [F,  sn(n) >=   sumtemp - lamda * norm(  LC(j,:) - His(n,:) ,1) ];
    end
    sumtemp = 0;
end
%F = [F,  lamda >=0,yij>=0,yij<=1];
F = [F, lamda >=0 ];
ops = sdpsettings('solver','cplex','verbose',2);
ops.cplex.benders.strategy = 3 ;
%     ops.gurobi.MIPGap = 0.005;
%     ops.gurobi.MIPGapAbs= 0.005;
%     ops.gurobi.OptimalityTol= 0.005;
solmp=solvesdp(F,Obj,ops );
solmp.info;
if isempty(strfind(solmp.info,'Successfully'))
    warning('不可行！更改最大theta的值！！')
end
Obj = value(Obj)
beta = value(lamda);
yij =  value(yij);
xi =  value(xi)
an =  value(sn) ;

toc

% 
figure(2)
open = find(xi ==1)  ;
plot( DC(open,1),DC(open,2) , 'rs','MarkerFaceColor','k')
hold on
plot(LC(:,1), LC(:,2), 'k.')
axis equal
hold on
plot(His(:,1), His(:,2), 'bo')
hold on
plot (DC(:,1), DC(:,2), 's')
% 
% for i = 1:size(DC,1)
%     for j = 1:size(LC,1)
%         if yij(i,j) == 1
%             plot([DC(i,1), LC(j,1)], [DC(i,2), LC(j,2)] ,'g-') ;
%             hold on
%         end
%     end
% end




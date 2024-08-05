clc;clear;
tic
load His ;
 % load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
%  [DC] = Distance(DC) ;
b = 10;  % The max number of opened DC
aa = 29;  % number of demand points in x-axes
bb = 29  ;% ...in y-axes
MultiFactor = 1.3 ;
aaa = 7 ;
bbb = 4 ;
[His] = Distance(His) ;
% minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
% maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
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

EmProba = 1/size(His,1) ;
% µÃµ½ Z_0
z_0 = getZ0(EmProba, DC, His, b) ;
tau = z_0 * MultiFactor ;

xi = binvar(size(DC,1),1,'full');
yij = binvar(size(DC,1),size(LC,1),'full'); %
sn = sdpvar(size(His,1),1,'full');
k = sdpvar(1,1,'full');

%% Objective
Obj = k;

%% Constraints
F = [] ;
F = [F, EmProba*sum(sn) <= tau ] ;

sumtemp = 0;
for j = 1:size(LC,1)
    for i = 1:size(DC,1)
        sumtemp = sumtemp + yij(i,j) * norm(DC(i,:) - LC(j,:),2) ;
    end
    for n = 1:size(His,1)
        F = [F,  sn(n) >= sumtemp - k * norm(  LC(j,:) - His(n,:) ,1) ];
    end
    sumtemp = 0;
end

F = [F, sum(xi) == b] ;

for j = 1:size(LC,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end
F = [F, k>=0, sn>=0] ;


ops = sdpsettings('solver','gurobi','verbose',3);
%ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,Obj,ops );
solmp.info;
Obj = value(Obj)   
xi =  value(xi)' 
sn =  value(sn) ;
yij = value(yij) ;
k = value(k) ;
%xx = xi'  ;

  [Obj_Value1] = Func_Test_RSM0001(x')
   [Obj_Value5] = Func_Test_RSM0005(x')

toc






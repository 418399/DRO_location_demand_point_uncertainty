  clc;clear; 
 
load Sichuan_DC;
%Test_loca = Distance(Test_loca);
b = 10;   % The max number of opened Test_loca
 
Test_loca = [101.85	32.24
104.82	29.59
104.79	29.55
104.77	28.43
104.89	28.37
104.86	28.2
105	32.27
103	30.3
102.9	30.1
99.4	31
105.3	32.6
100.9	31.3
105.7	30.3
105.6	32.7
105.3	32.5
103.4	31.2
104.5	31.7
103.5	31.3
103.9	31.5
104.1	31.3];
 
 
 

yij = binvar(size(Sichuan_DC,1),size(Test_loca,1),'full'); %
xi = binvar(size(Sichuan_DC,1),1,'full');

%% Objective
Obj = 0;
for j = 1:size(Test_loca,1)
    for i = 1:size(Sichuan_DC,1)
        Obj = Obj +  ( yij(i,j) * norm(Sichuan_DC(i,:) - Test_loca(j,:) ,2)  ) ;
    end
end
 

%% Constraints
F = [];
% number of opened facilities
F = [F, sum(xi) == b] ;

% allocation decision
for j = 1:size(Test_loca,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(Sichuan_DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end


ops = sdpsettings('solver','','verbose',2 );
solmp=solvesdp(F,Obj,ops );

solmp.info;
Obj = value(Obj)


yij =  value(yij)  ;
xi =  value(xi) ;
 [averagedis] = sichuantest(xi,Sichuan_DC) 













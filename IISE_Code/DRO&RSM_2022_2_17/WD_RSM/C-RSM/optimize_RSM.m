function [xi,yij, sn, k, Obj] = optimize_RSM(EmProba, DC, Sbar, His, b, tau)


xi = binvar(size(DC,1),1,'full');
yij = sdpvar(size(DC,1),size(Sbar,1),'full'); %
sn = sdpvar(size(His,1),1,'full');
k = sdpvar(1,1,'full');

%% Objective
 Obj = k;

 %% Constraints
 F = [] ;
 
 F = [F, EmProba*sum(sn) <= tau ] ;

 sumtemp = 0;
for j = 1:size(Sbar,1)
    for i = 1:size(DC,1)
        sumtemp = sumtemp + yij(i,j) * norm(DC(i,:) - Sbar(j,:),2) ;
    end
    for n = 1:size(His,1)
        F = [F,  sn(n) >=   sumtemp - k * norm(  Sbar(j,:) - His(n,:) ,1) ];
    end
    sumtemp = 0;
end

 F = [F, sum(xi) == b] ;

for j = 1:size(Sbar,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end
F = [F, k>=0, yij >=0,yij<=1] ;


ops = sdpsettings('solver','gurobi','verbose',0);
%ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,Obj,ops );
 solmp.info;
Obj = value(Obj)  ; %记录，松弛解为下界。
xi =  value(xi) ;
sn =  value(sn) ;
yij = value(yij) ;
k = value(k) ;



end


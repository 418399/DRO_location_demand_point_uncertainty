function [xi,an, beta,Obj  ] = optimize(theta,b,His,DC,Sbar,EmProba)

an = sdpvar(size(His,1),1,'full');
beta = sdpvar(1,1,'full');
xi = binvar(size(DC,1),1,'full');
yij = sdpvar(size(DC,1),size(Sbar,1),'full'); %
%% Objective
%Obj = 0;
Obj = beta * theta + EmProba*sum(an) ;
%% Constraints
F = [];
% number of opened facilities
F = [F, sum(xi) == b] ;
% allocation decision
for j = 1:size(Sbar,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end
% Dual constraints
sumtemp = 0;
for j = 1:size(Sbar,1)
    for i = 1:size(DC,1)
        sumtemp = sumtemp + yij(i,j) * norm(DC(i,:) - Sbar(j,:),2) ;
    end
    for n = 1:size(His,1)
        F = [F,  an(n) >=   sumtemp - beta * norm(  Sbar(j,:) - His(n,:) ,1) ];
    end
    sumtemp = 0;
 
end

%F = [F, beta >=0 ];
F = [F, beta >=0,yij>=0,yij<=1];
ops = sdpsettings('solver','cplex','verbose',3);
% ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,Obj,ops );
%solmp.info;
Obj = value(Obj)  ; %记录，松弛解为下界。
xi =  value(xi) ;
an =  value(an) ;
beta = value(beta) ;
yij = value(yij) ;

    if isempty(strfind(solmp.info,'Successfully'))
         warning('不可行！更改最大theta的值！！')
         warinngg = 1 ;
    end
end


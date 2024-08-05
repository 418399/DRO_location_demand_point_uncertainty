function [Obj, x] = getZ0_distance(EmProba, DC, His,b  )
 % His = Distance(His);
yij = binvar(size(DC,1),size(His,1),'full'); %
xi = binvar(size(DC,1),1,'full');

%% Objective
Obj = 0;
for j = 1:size(His,1)
    for i = 1:size(DC,1)
      %  Obj = Obj +  ( yij(i,j) * norm(DC(i,:) - His(j,:) ,2)  ) ;
      nor = distance(DC(i,1),DC(i,2),His(j,1),His(j,2))/180*pi*6371 ;
      Obj = Obj +  ( yij(i,j) * nor  ) ;
    end
end
Obj = EmProba * Obj ;  % Probability equality

%% Constraints
F = [];
% number of opened facilities
F = [F, sum(xi) == b] ;

% allocation decision
for j = 1:size(His,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end


ops = sdpsettings('solver','','verbose', 0);
solmp=solvesdp(F,Obj,ops );
 
Obj = value(Obj)  ;
x = value(xi) ;










end


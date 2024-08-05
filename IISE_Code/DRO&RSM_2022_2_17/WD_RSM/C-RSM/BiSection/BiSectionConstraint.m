function  [Status,xi] = BiSectionConstraint(maxb,b,His,LC,DC,tau,DiffertKK,K)


xi = binvar(size(DC,1),1,'full');
yij1 = binvar(size(DC,1),size(DiffertKK{1,2},1),'full'); %
yij2 = binvar(size(DC,1),size(DiffertKK{2,2},1),'full'); %
yij3 = binvar(size(DC,1),size(DiffertKK{3,2},1),'full'); %
% yij4 = binvar(size(DC,1),size(DiffertKK{4,2},1),'full'); %
% yij5 = binvar(size(DC,1),size(DiffertKK{5,2},1),'full'); %
% alpha1n = sdpvar(size(DiffertKK{1,1},1), 1,'full');
% alpha2n = sdpvar(size(DiffertKK{2,1},1), 1,'full');
% alpha3n = sdpvar(size(DiffertKK{3,1},1), 1,'full');
% alpha4n = sdpvar(size(DiffertKK{4,1},1), 1,'full');
% alpha5n = sdpvar(size(DiffertKK{5,1},1), 1,'full');


F = [];  % Initialize the costraints
F = [F, sum(xi) == b] ;
Obj1 =0;
Obj2 =0;
Obj3 =0;
% Obj4 =0;
% Obj5 =0;

for k = 1:K
    pk(k) = size(DiffertKK{k,1},1) / size(His,1) ;
    pro_nk(k) = 1/size(DiffertKK{k,1},1) ;  % the probability of N_k
    if k==1 % Constraints
        for n = 1:size(DiffertKK{k,1},1)  %His
            for j = 1:size(DiffertKK{k,2},1)   % LC
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij1(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij1(i,j) <= xi(i)  ];  % allocation decision
                end
                Obj1 = Obj1 + pk(k) * ( temp - maxb * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) )   ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij1(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==2 % Constraints
        for n = 1:size(DiffertKK{k,1},1)  %His
            for j = 1:size(DiffertKK{k,2},1)   % LC
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij2(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij2(i,j) <= xi(i)  ];  % allocation decision
                end
                Obj2 = Obj2 + pk(k) * ( temp - maxb * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) )   ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij2(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==3 % Constraints
        for n = 1:size(DiffertKK{k,1},1)  %His
            for j = 1:size(DiffertKK{k,2},1)   % LC
                temp = 0;
                for i = 1:size(DC,1)
                    temp = temp + yij3(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij3(i,j) <= xi(i)  ];  % allocation decision
                end
                Obj3 = Obj3 + pk(k) * ( temp - maxb * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) )  ;   % every j  , k,  n has a corresponding constraint
                F = [F, sum(yij3(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    %     if k==4 % Constraints
    %         for n = 1:size(DiffertKK{k,1},1)  %His
    %             for j = 1:size(DiffertKK{k,2},1)   % LC
    %                 temp = 0;
    %                 for i = 1:size(DC,1)
    %                     temp = temp + yij4(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij4(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 Obj4 = Obj4 + pk(k) * ( temp - maxb * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) )   ;   % every j  , k,  n has a corresponding constraint
    %                 F = [F, sum(yij4(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    %
    %
    %     if k==5 % Constraints
    %         for n = 1:size(DiffertKK{k,1},1)  %His
    %             for j = 1:size(DiffertKK{k,2},1)   % LC
    %                 temp = 0;
    %                 for i = 1:size(DC,1)
    %                     temp = temp + yij5(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij5(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 Obj5 = Obj5 + pk(k) * ( temp - maxb * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) )   ;   % every j  , k,  n has a corresponding constraint
    %                 F = [F, sum(yij5(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    
end

Obj = [Obj1 Obj2 Obj3 ] * pro_nk'    ;
ops = sdpsettings('solver','gurobi','verbose',0 );
%ops.cplex.benders.strategy = 3 ;
solmp=solvesdp(F,-Obj,ops );
%solmp.info;

xi = value(xi) ;
Obj = value(Obj) ;

if Obj <= tau
    Status =1;
else
    Status =0;
    
end






end

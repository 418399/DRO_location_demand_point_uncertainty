clc;clear;tic
load His ;
 
His = Sourcedata1(1:20,:);  % Historical data point
 His = Distance(His);
b = 10;   % The max number of opened DC
aa = 29;  % number of demand points in x-axes
bb = 29  ;% ...in y-axes
aaa = 7 ;
bbb = 4 ;
K=  3 ;  % @@@@@@每次改K后，记得改中心点的load文件，循环中的if k==？ 和 约束中的>=0!!!!!!!!!

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
a = [];  % Find the MAX theta
for j = 1:size(LC,1)
    for n = 1:size(His,1)
        a(j,n) = norm(His(n,:) - LC(j,:),1  ) ;
    end
end
a = sum(a,2) ;
maxtheta = max(a) * EmProba;
%  No_initialLC = round(0.2* size(LC,1) ) ;
%   random_Num = randperm(size(LC,1),No_initialLC) ;
%    [idx,Center] = kmeans(His,K,'Distance','sqeuclidean', 'MaxIter',10000 );  % K-means cluster

load idx;
load Center ;
His = [His idx] ;   % add the cluster number of each points in the 3th column
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
Vcells = Voronoi(Center ,area)  ;
hold on
plot(His(:,1),His(:,2),'v')
% axis equal
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



t = sdpvar(1,1) ;

xi = binvar(size(DC,1),1,'full');
yij1 = binvar(size(DC,1),size(DiffertKK{1,2},1),'full'); %
yij2 = binvar(size(DC,1),size(DiffertKK{2,2},1),'full'); %
yij3 = binvar(size(DC,1),size(DiffertKK{3,2},1),'full'); %
% yij4 = binvar(size(DC,1),size(DiffertKK{4,2},1),'full'); %
% yij5 = binvar(size(DC,1),size(DiffertKK{5,2},1),'full'); %
%     yij6 = binvar(size(DC,1),size(DiffertKK{6,2},1),'full'); %
%     yij7 = binvar(size(DC,1),size(DiffertKK{7,2},1),'full'); %

beta = sdpvar(1,1) ;
alpha1n = sdpvar(size(DiffertKK{1,1},1), 1,'full');
alpha2n = sdpvar(size(DiffertKK{2,1},1), 1,'full');
alpha3n = sdpvar(size(DiffertKK{3,1},1), 1,'full');
% alpha4n = sdpvar(size(DiffertKK{4,1},1), 1,'full');
% alpha5n = sdpvar(size(DiffertKK{5,1},1), 1,'full');
%     alpha6n = sdpvar(size(DiffertKK{6,1},1), 1,'full');
%     alpha7n = sdpvar(size(DiffertKK{7,1},1), 1,'full');

F = [];  % Initialize the costraints
F = [F, sum(xi) == b] ;

%% Obj

temp = 0;
for k = 1:K
    pk(k) = size(DiffertKK{k,1},1) / size(His,1) ;
    pro_nk(k) = 1/size(DiffertKK{k,1},1) ;  % the probability of N_k
    
    if k==1 % Obj
        for n = 1:size(DiffertKK{k,1},1)
            temp = temp + pk(k)/pro_nk(k) * alpha1n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                tempC = 0;
                for i = 1:size(DC,1)
                    tempC = tempC+ yij1(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij1(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha1n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
                F = [F, sum(yij1(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==2 % Obj
        for n = 1:size(DiffertKK{k,1},1)
            temp = temp + pk(k)/pro_nk(k) * alpha2n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                tempC = 0;
                for i = 1:size(DC,1)
                    tempC = tempC+ yij2(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij2(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha2n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
                F = [F, sum(yij2(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    if k==3 % Obj
        for n = 1:size(DiffertKK{k,1},1)
            temp = temp + pk(k)/pro_nk(k) * alpha3n(n) ;
            for j = 1:size(DiffertKK{k,2},1)
                tempC = 0;
                for i = 1:size(DC,1)
                    tempC = tempC+ yij3(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
                    F = [F, yij3(i,j) <= xi(i)  ];  % allocation decision
                end
                F = [F, alpha3n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
                F = [F, sum(yij3(:,j)) == 1] ;  % allocation decision
            end
        end
    end
    
    %     if k==4 % Obj
    %         for n = 1:size(DiffertKK{k,1},1)
    %             temp = temp + pk(k)/pro_nk(k) * alpha4n(n) ;
    %             for j = 1:size(DiffertKK{k,2},1)
    %                 tempC = 0;
    %                 for i = 1:size(DC,1)
    %                     tempC = tempC+ yij4(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij4(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 F = [F, alpha4n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
    %                 F = [F, sum(yij4(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    %
    %     if k==5 % Obj
    %         for n = 1:size(DiffertKK{k,1},1)
    %             temp = temp + pk(k)/pro_nk(k) * alpha5n(n) ;
    %             for j = 1:size(DiffertKK{k,2},1)
    %                 tempC = 0;
    %                 for i = 1:size(DC,1)
    %                     tempC = tempC+ yij5(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij5(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 F = [F, alpha5n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
    %                 F = [F, sum(yij5(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    %
    %     if k==6 % Obj
    %         for n = 1:size(DiffertKK{k,1},1)
    %             temp = temp + pk(k)/pro_nk(k) * alpha6n(n) ;
    %             for j = 1:size(DiffertKK{k,2},1)
    %                 tempC = 0;
    %                 for i = 1:size(DC,1)
    %                     tempC = tempC+ yij6(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij6(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 F = [F, alpha6n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
    %                 F = [F, sum(yij6(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    %
    %     if k==7 % Obj
    %         for n = 1:size(DiffertKK{k,1},1)
    %             temp = temp + pk(k)/pro_nk(k) * alpha7n(n) ;
    %             for j = 1:size(DiffertKK{k,2},1)
    %                 tempC = 0;
    %                 for i = 1:size(DC,1)
    %                     tempC = tempC+ yij7(i,j) * norm( DC(i,:) - DiffertKK{k,2}(j,:) ,2) ;% calculate the distance between DC and demand points
    %                     F = [F, yij7(i,j) <= xi(i)  ];  % allocation decision
    %                 end
    %                 F = [F, alpha7n(n) + t >=  tempC - beta * norm( DiffertKK{k,2}(j,:) -  DiffertKK{k,1}(n,:) ,1) ]  ;
    %                 F = [F, sum(yij7(:,j)) == 1] ;  % allocation decision
    %             end
    %         end
    %     end
    
    k
    
end
F = [F, beta>=0, t>=0, alpha1n>=0,alpha2n>=0, alpha3n>=0 ]  ;
save_res = {};
%  ,alpha4n>=0,alpha5n>=0
% etaS = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
% multiplierS = [0.02 0.06 0.1 0.14 0.18	0.2];  %The factor to control the parameter \theta
etaS =  [0.1  ];
multiplierS = [0.02 0.06 0.1 0.14 0.18];  %The factor to control the parameter \theta
for etakkk = 1:size(etaS,2)
    eta = etaS(etakkk) ;
    for eee = 1:size(multiplierS,2)
        theta =multiplierS(eee) * maxtheta ;
        Obj = 0;
        Obj = Obj + t + (1/eta) * theta * beta ;
        Obj = Obj + (1/eta) * temp ;
        
        ops = sdpsettings('solver','cplex','verbose',0 );
        solmp=solvesdp(F,Obj,ops );
        solmp.info
        x = value(xi);
        ts = value(t) ;
        SaObj = value(Obj) ;
        abc = toc;
        [averagedis] = CVaR_DRO_Test(x,DC,His);
        save_res(etakkk,eee) = {[x; ts; averagedis; abc] };
        
        
    end
end


















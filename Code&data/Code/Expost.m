clc;clear;

load His ;
load Test_points_35;
His = Sourcedata1(1:20,:);  % Historical data point
[His] = Distance(His) ;
aaa = 7 ;
bbb = 4 ;
b = 10;  % The max number of opened DC

minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
% minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
% maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
lengthY = maxlo(2) - minlo(2) ;
 
%% DC
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end

Test_poi = [];
for ii = 1:10
    Test_poi = [Test_poi;  Test_points(:,2*ii-1:2*ii) ] ;
end



yij = binvar(size(DC,1),size(Test_poi,1),'full'); %
xi = binvar(size(DC,1),1,'full');
%% Objective
Obj = 0;
for j = 1:size(Test_poi,1)
    for i = 1:size(DC,1)
        Obj = Obj +  ( yij(i,j) * norm(DC(i,:) - Test_poi(j,:) ,2)  ) ;
    end
end
 
%% Constraints
F = [];
F = [F, sum(xi) == b] ; % number of opened facilities

% allocation decision
for j = 1:size(Test_poi,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end

ops = sdpsettings('solver','','verbose',0 );
solmp=solvesdp(F,Obj,ops );
solmp.info;
y =  value(yij);
x =  value(xi)  ;
%Objective = value(Obj)  ;



for i=1:size(DC,1)  % DC到所有DP的2范数
    for j = 1:size(Test_poi,1)
        save_2_norm(i,j) = norm( DC(i,:) -Test_poi(j,:) ,2)  ;
    end
end
temp_save_2_norm = save_2_norm( find(x==1),:)  ;
[min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % 对于每个z，找到与其最近的DC及距离
Save_results = min_value ;
Save_CVaR_Result = CVaRR(Test_poi, DC, x)  ;

% Print results
MAX_test = max(Save_results) ;
MEAN_test = mean(Save_results);
Save_CVaR_Result;
Print_results = [MAX_test ;  MEAN_test ; Save_CVaR_Result] 




 






% open = find(xi ==1)  ;
% plot( DC(open,1),DC(open,2) , 'ro','MarkerFaceColor','r')
% hold on
% %plot(LC(:,1), LC(:,2), 'kv')
% hold on
% plot(Test_loca(:,1), Test_loca(:,2), 'bo')
% hold on
% plot (DC(:,1), DC(:,2), 's')
% axis equal



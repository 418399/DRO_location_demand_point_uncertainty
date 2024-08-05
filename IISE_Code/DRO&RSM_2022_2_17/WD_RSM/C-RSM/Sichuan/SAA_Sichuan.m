  clc;clear; 
load Sichuan_data;
load Sichuan_DC;
%Sichuan_data = Distance(Sichuan_data);
b = 10;   % The max number of opened Sichuan_data
% aa = 39;  % number of demand points in x-axes
% bb = 39  ;% ...in y-axes
% K=  7 ;
EmProba = 1/size(Sichuan_data,1) ;
%  minlo = [min([Sichuan_data(:,1) ]) min([Sichuan_data(:,2) ])] ;
% maxlo = [max([Sichuan_data(:,1) ]) max([Sichuan_data(:,2)])];
% % minlo = [min([Sichuan_data(:,1) ; Sichuan_Sichuan_data(:,1)]) min([Sichuan_data(:,2) ;Sichuan_Sichuan_data(:,2)])] ;
% % maxlo = [max([Sichuan_data(:,1) ;Sichuan_Sichuan_data(:,1)]) max([Sichuan_data(:,2);Sichuan_Sichuan_data(:,2)])];
% lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
% lengthY = maxlo(2) - minlo(2) ;

% Disc_X = [minlo(1):lengthX/aa:maxlo(1)]  ;  % get the elements of the x and y  axes.
% Disc_Y = [minlo(2):lengthY/bb:maxlo(2)] ;

% LC = [];   %  Get the square lattice  locations
% for i = 1:size(Disc_X,2)
%     for j = 1:size(Disc_Y,2)
%         LC = [LC ; Disc_X(i) , Disc_Y(j) ] ;
%     end
% end

yij = binvar(size(Sichuan_DC,1),size(Sichuan_data,1),'full'); %
xi = binvar(size(Sichuan_DC,1),1,'full');

%% Objective
Obj = 0;
for j = 1:size(Sichuan_data,1)
    for i = 1:size(Sichuan_DC,1)
        Obj = Obj +  ( yij(i,j) * norm(Sichuan_DC(i,:) - Sichuan_data(j,:) ,2)  ) ;
    end
end
Obj = EmProba * Obj ;  % Probability equality

%% Constraints
F = [];
% number of opened facilities
F = [F, sum(xi) == b] ;

% allocation decision
for j = 1:size(Sichuan_data,1)
    F = [F, sum( yij(:,j) ) ==1 ];
    for i = 1:size(Sichuan_DC,1)
        F = [F, yij(i,j) <= xi(i)  ];
    end
end


ops = sdpsettings('solver','','verbose',2 );
solmp=solvesdp(F,Obj,ops );

solmp.info;
Obj = value(Obj)


yij =  value(yij);
xi =  value(xi)  ;

 [averagedis] = sichuantest(xi,Sichuan_DC) 

% open = find(xi ==1)  ;
% plot( Sichuan_data(open,1),Sichuan_data(open,2) , 'ro','MarkerFaceColor','r')
% hold on
% %plot(LC(:,1), LC(:,2), 'kv')
% hold on
% plot(Sichuan_data(:,1), Sichuan_data(:,2), 'bo')
% hold on
% plot (Sichuan_data(:,1), Sichuan_data(:,2), 's')
% axis equal

 

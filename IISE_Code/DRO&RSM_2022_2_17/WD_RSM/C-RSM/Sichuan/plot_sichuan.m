clc;clear;

load Sichuan_data;
load Sichuan_DC;
tempSichuan_DC = Sichuan_DC;
tempSichuan_DC(10,:) = [];
tempSichuan_DC(15,:) = [];
tempSichuan_DC(14,:) = [];
tempSichuan_DC(36,:) = [];
minlo = [min([Sichuan_data(:,2); tempSichuan_DC(:,2) ]) min([Sichuan_data(:,1); tempSichuan_DC(:,1)    ])] ;
maxlo = [max([Sichuan_data(:,2) ;tempSichuan_DC(:,2) ]) max([Sichuan_data(:,1) ;tempSichuan_DC(:,1)   ])];
load idxsichuan7;
load Centersichuan7 ;

tSichuan_data = [Sichuan_data(:,2), Sichuan_data(:,1)]  ;
Sichuan_data = [tSichuan_data idx ] ;

%area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
tarea = [97.201667	25.9; 97.201667	34.2;108.718375	34.2;108.718375	25.9] ;

area = [tarea(:,2), tarea(:,1)] ;

Vcells = Voronoiplot(Center ,area)  ;
 
hold on
 demand =   plot(Sichuan_data(:,1),Sichuan_data(:,2),'r.', 'MarkerSize',10)
x = [0	1	1	0	1	0	0	0	0	0	1	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0
];
 for i  =1:size(x,2)
    if x(1,i)==1
     sichuandcopen = plot(Sichuan_DC(i,2),Sichuan_DC(i,1),'ko','MarkerFaceColor','k','MarkerSize',5);
    end
 end

  for i  =1:size(tempSichuan_DC,1)
      sichuandc = plot(tempSichuan_DC(i,2),tempSichuan_DC(i,1),'ko','MarkerSize',5);
  end
 
 trainingdata = [105.6	32.7
105.50 	32.70 
105.3	32.6
105.40 	32.60 
105.40 	32.60 
105.3	32.5
105.00 	32.30 
104.90 	32.30 
104.80 	32.30 
105	32.27
101.85	32.24
104.70 	32.20 
104.5	31.7
104.00 	31.60 
104.20 	31.60 
104.00 	31.50 
103.9	31.5
103.80 	31.40 
103.90 	31.30 
103.30 	31.30 
103.50 	31.30 
103.60 	31.30 
103.80 	31.30 
104.1	31.3
103.50 	31.20 
103.4	31.2
103.70 	31.10 
99.4	31
102.90 	30.30 
102.9	30.1
103.30 	31.00 
104.82	29.59
104.79	29.55
104.77 	28.43 
104.89	28.37
105.03	28.22
104.86	28.2
99.30 	28.20 
];

 
 hold on
 train = plot(trainingdata(:,2),trainingdata(:,1),'ks');
 


axis equal
% legend([train demand sichuandc sichuandcopen ], 'Training points' , 'Test points', 'Candidate DC', 'Opened DC'  );
     xlabel('longitude')
ylabel('latitude')
ylim([97.3,107])
xlim([26,34])
grid on


%% 设施开放数量对平均距离的影响

DIS = [38.3186	30.58915777	25.50916828	22.0351691	21.39239269];
x = [10 20 30 40 50] ;
plot(x,DIS,'k-');
 xlabel('Number of opened facilities')
ylabel('Average distance (km)')







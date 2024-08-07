
clc;clear;
load His ;
load DC;
%   load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
Expost = [0	0	0	0	1	0	1	0	0	1	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	1	0	0
    0	0	0	0	1	0	0	0	1	1	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	0
    0	0	0	1	1	0	0	1	0	0	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	1	0	0
    1	0	0	1	1	0	0	1	0	0	0	1	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	1	0	0
    ];
DRO = [0	0	0	0	1	0	1	0	1	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	1	0	0
    0	0	0	0	1	0	1	0	1	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	1	0	0
    0	0	0	1	1	0	1	0	0	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	1	0	0
    0	0	0	1	1	0	1	0	0	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	1	0	0
    ];

RS = [0	0	0	0	1	0	1	0	0	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	1	1	0
    0	0	0	0	1	0	1	0	1	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	0
    0	0	0	1	1	0	1	0	0	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	1	0	0
    0	0	0	1	1	0	1	0	0	1	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	1	0	0
    ];



figure () ;
hold on;
for i  =1:size(Expost,2)
    if Expost(1,i)==1
        plot(DC(i,1),DC(i,2),'ko','MarkerSize',10);
    end
    if DRO(1,i)==1
        plot(DC(i,1),DC(i,2),'kv');
    end
    if RS(1,i)==1
        plot(DC(i,1),DC(i,2),'k.');
    end
end
axis equal

figure () ;
hold on;
for i  =1:size(Expost,2)
    if Expost(2,i)==1
        plot(DC(i,1),DC(i,2),'ko','MarkerSize',10);
    end
    if DRO(2,i)==1
        plot(DC(i,1),DC(i,2),'kv');
    end
    if RS(2,i)==1
        plot(DC(i,1),DC(i,2),'k.');
    end
end
axis equal

figure () ;
hold on;
for i  =1:size(Expost,2)
    if Expost(3,i)==1
        plot(DC(i,1),DC(i,2),'ko','MarkerSize',10);
    end
    if DRO(3,i)==1
        plot(DC(i,1),DC(i,2),'kv');
    end
    if RS(3,i)==1
        plot(DC(i,1),DC(i,2),'k.');
    end
end
axis equal

figure () ;
hold on;
for i  =1:size(Expost,2)
    if Expost(4,i)==1
        plot(DC(i,1),DC(i,2),'ko','MarkerSize',10);
    end
    if DRO(4,i)==1
        plot(DC(i,1),DC(i,2),'kv');
    end
    if RS(4,i)==1
        plot(DC(i,1),DC(i,2),'k.');
    end
end
axis equal


%% 不同参数取值下，不同BW的表现
DRO = [1.4725	1.9963	2.5682	3.05	3.51	3.6891
1.4400 	1.9398 	2.4577 	2.9515 	3.4067 	3.6586 
1.4400 	1.8398 	2.4277 	2.9315 	3.4067 	3.6586 
1.4400 	1.8498 	2.4177 	2.9015 	3.4067 	3.6586 
1.4400 	1.8770 	2.4577 	2.8287 	3.1500 	3.6586 
1.5338 	1.9664 	2.4819 	2.9315 	3.2867 	3.4902 
1.6072 	2.0197 	2.5655 	3.0200 	3.365625925	3.549626682
1.8468 	2.6008 	2.972516902	3.218652623	3.386509042	3.560208927
1.8947 	2.90628594	3.008795471	3.236714559	3.40739216	3.570791172
2.193163441	2.943191893	3.045074041	3.254776495	3.428275277	3.581373418
2.429923785	2.980097845	3.08135261	3.272838431	3.449158395	3.591955663
2.904642013	3.017003797	3.11763118	3.290900367	3.470041513	3.602537908
2.93287562	3.05390975	3.15390975	3.308962303	3.49092463	3.613120154
];

x = [0 0.02	0.04	0.06	0.08	0.1	0.12	0.14	0.16	0.18	0.2	0.3	0.4];

figure ()
hold on
bw5 = plot(x,DRO(:,1), 'kd-')
bw15 = plot(x,DRO(:,2), 'kv-')
bw25 = plot(x,DRO(:,3), 'ks-')
bw35 = plot(x,DRO(:,4), 'ko-')
bw45 = plot(x,DRO(:,5), 'k^-')
bw55 = plot(x,DRO(:,6), 'k*-')
 legend([bw5 bw15 bw25 bw35 bw45 bw55], 'B_w=0.5', 'B_w=1.5', 'B_w=2.5', 'B_w=3.5', 'B_w=4.5', 'B_w=5.5' );

 
 
 RS = [1.4725	1.9963	2.5682	3.05	3.51	3.6891
1.4598 	1.9052 	2.5297 	3.0261 	3.4045 	3.6171 
1.5115 	1.8613 	2.4100 	2.8936 	3.2912 	3.5103 
1.5127 	1.8833 	2.2904 	2.8762 	3.2679 	3.4558 
1.5127 	1.8996 	2.4114 	2.7610 	3.2176 	3.4558 
1.5127 	1.8996 	2.4114 	2.8287 	3.1300 	3.4036 
1.5338 	1.9664 	2.4219 	2.8762 	3.2076 	3.3602 
1.6072 	2.0197 	2.5655 	3.0200 	3.3368 	3.5658 
1.6806 	2.0730 	2.7091 	3.2112 	3.5436 	3.7915 
];
 RSx = [1	1.05	1.1	1.2	1.3	1.4	1.5	1.6	1.7];
 figure ()
hold on
bw5 = plot(RSx,RS(:,1), 'kd-')
bw15 = plot(RSx,RS(:,2), 'kv-')
bw25 = plot(RSx,RS(:,3), 'ks-')
bw35 = plot(RSx,RS(:,4), 'ko-')
bw45 = plot(RSx,RS(:,5), 'k^-')
bw55 = plot(RSx,RS(:,6), 'k*-')
 legend([bw5 bw15 bw25 bw35 bw45 bw55], 'B_w=0.5', 'B_w=1.5', 'B_w=2.5', 'B_w=3.5', 'B_w=4.5', 'B_w=5.5' );

 
 
 
 
 

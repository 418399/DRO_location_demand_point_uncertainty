
% 离散与连续求解时间及精度对比
jingdu =  [10195.31514	10225.80998	10256.79097	10278.79097	10280.79097	10282.79097	10284.45573;
10339.86713	10339.86713	10339.86713	10339.86713	10339.86713	10339.86713	10339.86713 ]/1000;
 xxx = [1600 2500 3600 4900 6400 8100 10000]; 
plot(xxx,jingdu(1,:), '-') ;
hold on
plot(xxx,jingdu(2,:), '-') ;
% set(gca,'XTick',1:7);
%set(gca,'XTickLabel',{'40*40','50*50','60*60','70*70','80*80','90*90','100*100'})
xlabel('Scenario number')
ylabel('Distance')

 %% %
 figure()
time = [49.404323	71.538764	108.138259	124.182752	132.934106	207.956183	276.082012;
    2023 2023 2023 2023 2023 2023 2023];
 xxx = [1600 2500 3600 4900 6400 8100 10000]; 
plot(xxx,time(1,:), '-') ;
hold on
plot(xxx,time(2,:), '-') ;
% set(gca,'XTick',1:7);
% set(gca,'XTickLabel',{'40*40','50*50','60*60','70*70','80*80','90*90','100*100'})
xlabel('Scenario number')
ylabel('Solution time (s)')


%% 不同theta下，B-DRO与C-DRO的表现
  y01 = [5.88E+02 6.02E+02	6.29E+02	6.83E+02	7.82E+02 8.59E+02 9.59E+02 1.03E+03];
y02 = [9.52E+02 9.2778E+02	9.0555E+02	8.6296E+02	9.07E+02 9.25E+02 9.47E+02  9.68E+02]; 
  y01 =   y01/250;
  y02 = y02/250

x = [0.5 1 2 3 5 7 9 10] ;
plo01 = plot(x,y01, '--ks')
hold on
plo02 = plot(x,y02, '-ko')
hold on 
% plo3 = plot(0,5.203261492,'ko') ;
% plo4 = plot(0,5.8219,  'v')  ;
 legend([plo01 plo02  ], '\mu = 0.1', '\mu = 0.2' );
 %   legend([plo01 plo02 plo3 plo4], '\mu = 0.1', '\mu = 0.2', 'In-sample (\mu=0.1)',  'In-sample (\mu =0.2)');
xlabel('Bandwidth (km)')
ylabel('Average distance (km)')




%% 相同tau下，B-RS C-RS 不同测试案例下的表现:
brs = [1.50E+05	1.74E+05	1.86E+05	2.17E+05	2.52E+05];
crs = [1.50E+05	1.72E+05	1.82E+05	2.15E+05	2.51E+05] ;
brs = brs/100000;
crs = crs/100000;
x =[0.500 1  1.500 2 2.500 ];
p_b = plot(x,brs, '--ks')
hold on
p_c = plot(x,crs, '-ko')
grid on
 legend([p_b p_c  ], 'B-RS', 'C-RS' );
 xlabel('Bandwidth (km)')
ylabel('Average distance (km)')

% B_w较大时，2500，不同的tau下的表现
BRS = [2.52E+05	2.52E+05	2.42E+05	2.52E+05]; 
CRS = [2.51E+05	2.52E+05	2.42E+05	2.52E+05];
BRS = BRS/100000;
CRS = CRS/100000;
x = [1.1 1.3 1.5 1.7] ;
brs = plot(x,BRS, '--ks');
hold on
crs = plot(x,CRS,'-ko') ;
grid on
 legend([brs crs  ], 'B-RS', 'C-RS' );
 xlabel('$\tau /Z_0$')
ylabel('Average distance (km)')






%% \theta and tau 不同取值下，DRO 与 RSM 测试案例对比
CDRO_TEST = [1.46E+05	1.68E+05	1.85E+05	2.12E+05	2.64E+05	3.46E+05
1.51E+05	1.79E+05	1.96E+05	2.13E+05	2.42E+05	3.25E+05
1.63E+05	1.87E+05	2.15E+05	2.29E+05	2.69E+05	3.47E+05
1.87E+05	2.15E+05	2.22E+05	2.44E+05	2.71E+05	3.60E+05;
2.25E+05	2.32E+05		2.41E+05	2.53E+05	2.69E+05	3.52E+05
2.40E+05	2.68E+05	2.53E+05	2.81E+05	2.86E+05	3.59E+05];
CDRO_TEST = CDRO_TEST/100000 ;

CRS_TEST = [1.47E+05	1.72E+05	1.93E+05	2.22E+05	2.62E+05	3.65E+05	2.27E+05
                    1.50E+05	1.72E+05	1.82E+05	2.15E+05	2.51E+05	3.50E+05	2.20E+05
                    1.50E+05	1.74E+05	1.86E+05	2.13E+05	2.52E+05	3.43E+05	2.20E+05
                    1.51E+05	1.79E+05	1.96E+05	2.13E+05	2.42E+05	3.25E+05	2.18E+05
                    1.51E+05	1.79E+05	1.97E+05	2.19E+05	2.52E+05	3.15E+05	2.19E+05
                    1.87E+05	2.15E+05	2.22E+05	2.44E+05	2.71E+05	3.60E+05	2.50E+05];
CRS_TEST = CRS_TEST/100000;
x1 = [0.08	0.1	0.12	0.14 0.18	0.2];
x2 = [1.05 1.1	1.3	1.5	1.7	1.8] ;

figure()
hold on
cdro1000 = plot(x1,CDRO_TEST(:,2)', '-ks')
cdro2000 = plot(x1,CDRO_TEST(:,4)', '-kv')
%cdro5000 = plot(x1,CDRO_TEST(:,6)', '-ko')

%   legend([cdro1000 cdro2000  cdro5000], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
  legend([cdro1000 cdro2000  ], 'B_w=1km', 'B_w=2km' );
 xlabel('\mu')
ylabel('Average distance (km)')
hold off

figure()
hold on
crs1000 = plot(x2,CRS_TEST(:,2)', '-ks')
crs2000 = plot(x2,CRS_TEST(:,4)', '-kv')
crs_average = plot(x2,CRS_TEST(:,end), '--')
%crs5000 = plot(x2,CRS_TEST(:,6)', '-ko')
 legend([crs1000 crs2000 crs_average], 'B_w=1km', 'B_w=2km', 'Average');
%  legend([crs1000 crs2000 crs5000 ], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
 xlabel('\tau/Z_0')
ylabel('Average distance (km)')
hold off


%% 不同K下RSM 在Bw=5km下的表现，与SAA对比

RS=[3.47E+05	3.47E+05	3.25E+05	3.47E+05	3.47E+05 ]/100000; 
SAA = [3.59E+05	3.59E+05	3.59E+05	3.59E+05	3.59E+05]/100000 ;
DRO=[3.47E+05	3.47E+05	3.26E+05	3.47E+05	3.47E+05 ]/100000; 
xk = [1 2 3 4 5] ;

figure()
hold on
crs = plot(xk,RS, '-ks')
dro = plot(xk,DRO, '-ko')
saa = plot(xk,SAA, '-')
%crs5000 = plot(x2,CRS_TEST(:,6)', '-ko')
 legend([saa dro crs  ], 'SAA', 'C-DRO', 'C-RS' );
%  legend([crs1000 crs2000 crs5000 ], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
 xlabel('K (# Subregions)')
ylabel('Average distance (km)')
hold off



rs_saa_dro = [316697.5	316697.5	291137.5	316697.5	316697.5
322732.5	322732.5	322732.5	322732.5	322732.5
316697.5	316697.5	291137.5	303142.5	316697.5
]/100000 ;
xk = [1 2 3 4 5] ;

figure()
hold on
crs = plot(xk,rs_saa_dro(1,:), '-ks')
dro = plot(xk,rs_saa_dro(3,:), '-ko')
saa = plot(xk,rs_saa_dro(2,:), '-')
%crs5000 = plot(x2,CRS_TEST(:,6)', '-ko')
 legend([saa dro crs  ], 'SAA', 'C-DRO', 'C-RS' );
%  legend([crs1000 crs2000 crs5000 ], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
 xlabel('K (# Subregions)')
ylabel('Average distance (km)')


%%开放设施数量的测试表现
per_x_dro_rs = [6	5.208618265	5.08765905
8	4.634996004	4.362573356
10	3.849592517	3.849592517
12	3.689897156	3.689897156
14	3.645614912	3.472906345
16	3.317088293	3.307697188
];
figure()
hold on
crs = plot(per_x_dro_rs(:,1),per_x_dro_rs(:,3), '-ks')
dro = plot(per_x_dro_rs(:,1),per_x_dro_rs(:,2), '-ko')
 
%crs5000 = plot(x2,CRS_TEST(:,6)', '-ko')
 legend([ dro crs  ],  'C-DRO', 'C-RS' );
%  legend([crs1000 crs2000 crs5000 ], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
 xlabel('K (# Subregions)')
ylabel('Average distance (km)')



%% 不同问题规模下的求解时间：算法VS直接求解
figure()
hold on
timeA = [34.696556	94.820476	196.866414	203.533217	376.744813	597.48237] ;
timeDirec = [461 1400	3028	9000	19900	35920];
x = [900	1600	2500	3600	4900	6400];

al = plot(x,timeA, '-ks') ;
direcc = plot(x,timeDirec, '--ko')
 legend([al direcc],  'RCG algorithm', 'Direct solution method');
%  legend([crs1000 crs2000 crs5000 ], 'B_w=1km', 'B_w=2km', 'B_w=5km' );
 xlabel('Scenario number')
ylabel('Solution time (s)')



%% RS  DRO 在不同参数取值  以及 不同测试案例下的表现对比

DRO500 = [1.48E+05 1.48E+05 1.46E+05 1.46E+05 1.46E+05 1.51E+05 1.63E+05 1.87E+05  ]/100000;
DRO2500 = [2.68E+05 2.67E+05 2.64E+05 2.64E+05 2.64E+05 2.42E+05 2.69E+05 2.71E+05 ]/100000;

% RS500 = [1.47E+05	1.50E+05	1.50E+05	1.51E+05	1.51E+05	1.87E+05	  ]/100000;
% RS2500 = [2.62E+05	2.52E+05	2.51E+05	2.42E+05	2.52E+05	2.71E+05	 ]/100000;
RS500 = [1.48E+05 1.47E+05	1.50E+05	1.50E+05	1.50E+05	1.50E+05	1.51E+05	162820]/100000;
RS2500 = [2.68E+05 2.62E+05	2.51E+05	2.52E+05	2.52E+05	2.52E+05	2.42E+05	268530]/100000;

x = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14]; 
x2 = [1 1.05 1.1 1.2 1.3 1.4 1.5 1.6 ] ;


figure()
hold on
dro5_5_aver = plot(x,DRO500, '-ko') ;
dro5_25_aver = plot(x, DRO2500,'-ks');
% 画出mu=0时的测试解，作为比较。并画出expost解，作为benchmark，看做最佳解。
l5 = plot(x,[1.4504 1.4504 1.4504 1.4504 1.4504 1.4504 1.4504 1.4504],'k--')
l25 = plot(x,[2.2168 2.2168  2.2168 2.2168 2.2168 2.2168 2.2168 2.2168], 'k-')
 l55 = plot(0,1.48,'s', 'MarkerFaceColor','k')
l225 =  plot(0,2.68,'s', 'MarkerFaceColor','k')
  legend([dro5_5_aver dro5_25_aver l5 l25 l55], 'B_w=0.5km', 'B_w=2.5km', 'Ex post (B_w=0.5km)', 'Expost (B_w=2.5km)','Empirical');
 
xlabel('\mu')
ylabel('Average distance (km)')
grid on



figure()
hold on
 RS5_5_aver = plot(x2,RS500, '-ko' ) ;
 RS5_25_aver = plot(x2, RS2500,'-ks')
l5 = plot(x2,[1.4504 1.4504 1.4504 1.4504 1.4504 1.4504 1.4504 1.4504],'k-')
l25 = plot(x2,[2.2168  2.2168 2.2168 2.2168 2.2168 2.2168 2.2168 2.2168], 'k--')
l55=  plot(1,1.48,'s', 'MarkerFaceColor','k')
l225 =  plot(1,2.68,'s', 'MarkerFaceColor','k')
% 画出tau=1时的测试解，作为比较。并画出expost解，作为benchmark，看做最佳解。
  
   legend([RS5_5_aver RS5_25_aver l5 l25 l55], 'B_w=0.5km', 'B_w=2.5km', 'Ex post (B_w=0.5km)', 'Expost (B_w=2.5km)','Empirical');
 
 xlabel('\tau/Z_0')
ylabel('Average distance (km)')
grid on
xlim([1,1.6])


%% 不同B_w 下， SAA C-RO  C-DRO 的测试表现

resul_SAA_RS_DRO = [1.45E+05	1.71E+05	1.86E+05	2.17E+05	2.67E+05	3.59E+05	4.49E+05
1.50E+05	1.74E+05	1.86E+05	2.13E+05	2.52E+05	3.43E+05	404640
1.46E+05	1.68E+05	1.85E+05	2.12E+05	2.64E+05	3.46E+05	4.42E+05
1.45E+05	1.66E+05	1.74E+05	2.06E+05	2.22E+05	2.94E+05	3.23E+05]/100000;
X = [500	1000	1500	2000	2500	5000	10000]/1000;
figure()
hold on
 plot(X,resul_SAA_RS_DRO(1,:), '-', X, resul_SAA_RS_DRO(2,:),'-ks', X,resul_SAA_RS_DRO(3,:), '-kv',X,resul_SAA_RS_DRO(4,:), '--') ;
 
% 画出tau=1时的测试解，作为比较。并画出expost解，作为benchmark，看做最佳解。
% legend([al direcc],  'RCG algorithm', 'Direct solution method');
  legend(  'SAA', 'C-RS (\tau/Z_0=1.3)', 'C-DRO (\mu=0.02)' ,'Ex post');
 xlabel('B_w')
ylabel('Average distance (km)')


%%求解时间与初始子集大小的关系
time = [4928.49642	4000.612292	1721.705813	3808.581024	5828.49642	7959.200764	9988.149452
];
para = [0.03	0.06	0.1	0.2	0.3	0.4	0.5 ];
plot(para,time,'-ko');
 grid on
 xlabel('\omega ')
ylabel('Solution time (s)')


%% C-RS 在不同tau下，不同测试案例 B_w下 的测试表现

% crs11 = [1.47E+05	1.72E+05	1.93E+05	2.22E+05	2.62E+05	3.65E+05]/100000;
% crs15 = [1.51E+05 1.79E+05	1.96E+05	2.13E+05	2.42E+05	3.25E+05 ]/100000 ;
%  
% x =[0.5  1  1.5  2 2.5 5 ];
figure()
crs11 = [1.47E+05	1.72E+05	1.93E+05	2.22E+05	2.62E+05	 ]/100000;
crs15 = [1.51E+05 1.79E+05	1.96E+05	2.13E+05	2.42E+05	  ]/100000 ;
 
x =[0.5  1  1.5  2 2.5   ];
p_b = plot(x,crs11, '--ko')
hold on
p_c = plot(x,crs15, '-ks')
grid on
 legend([p_b p_c  ], '\tau=1.05', '\tau=1.5' );
 xlabel('Bandwidth (km)')
ylabel('Average distance (km)')


%% C-DRO 在不同tau下，不同测试案例 B_w下 的测试表现
% crs01 = [2.6816	2.8775	3.0898	3.0969	3.1507	4.3115];
% crs02 = [2.8461	2.8174	2.8293	2.9205	2.9511	3.7684];
%  
% x =[0.5  1  1.5  2 2.5 5 ];
figure()
crs001 = [1.46E+05	1.68E+05	1.85E+05	2.12E+05	2.64E+05 ]/100000;
crs01 = [1.51E+05	1.79E+05	1.96E+05	2.13E+05	2.42E+05]/100000;
 
x =[0.5  1  1.5  2 2.5   ];
p_b = plot(x,crs001, '--ko')
hold on
p_c = plot(x,crs01, '-ks')
grid on
 legend([p_b p_c  ], '\mu = 0.1', '\mu = 0.2' );
 xlabel('Bandwidth (km)')
ylabel('Average distance (km)')


%% 离散误差 
%   matix = [4854.1	4854.1	4854.1	4854.1	4854.1	4854.1 4854.1;
% 4718.4	4740	4797.6	4803.5	4810.4	4809 4836.4;
% 5213.37	5164.26	5151.15	5086.34	5022.53	4950.42 4907.11]/1000;
matix = [ 	4854.1	4854.1	4854.1	4854.1	4854.1	4854.1	4854.1;
 4740	4797.6	4803.5	4810.4	4.81E+03	4.84E+03	4845.9;
 5164.26	5151.15	5086.34	5022.53	4950.42	4907.11	4874.184]/1000;


x = [ 0.300	0.250	0.200	0.150	0.100 0.05 0.0];

figure()
hold on
a = plot(x,matix(1,:), 'k--')
b = plot(x,matix(2,:), 'ks-')
c = plot(x,matix(3,:), 'ko-')
grid on
 legend([a b c ], 'V_C (Continuous)', 'V_D (Discrete)', 'V_D+L(\delta)' );
 xlabel('L(\delta), (Scenario No. *10^2)')
ylabel('Distance (km)')
set(gca,'xticklabel',{  '0.02, (4000)'  '0.05, (640)' '0.1, (100)' '0.2, (49)'	'0.15, (64)' '0.25, (36)' '0.3, (25)'  })
 

 aa = [1.4400 	1.8770 	2.4577 	2.8287 	3.1900 	3.5586 
1.5127 	1.8596 	2.4114 	2.8000 	3.1300 	3.4036 
1.6531 	1.9732 	2.6033 	2.9754 	3.3111 	3.6666 
1.3111 	1.7200 	2.2133 	2.6210 	2.9899 	3.2421 ];
xx = [0.5	1.5	2.5	3.5	4.5	5.5];
figure ()
hold on
a = plot(xx,aa(1,:), 'ko-')
bb = plot(xx,aa(2,:), 'ks-')
cc = plot(xx,aa(3,:), 'k-')
dd = plot(xx,aa(4,:), 'k--')
grid on
 legend([a bb cc dd ], 'C-DRO (\theta=0.08\theta_{max})', 'C-RS (\tau/Z_0=1.4)', 'Empirical', 'Ex post'  );
 xlabel('B_w (km)')
ylabel('Minimized objective function value')




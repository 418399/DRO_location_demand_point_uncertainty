%% RS  DRO 在不同参数取值  以及 不同测试案例下的表现对比

DRO500 = [1.4725	1.4400 	1.4400 	1.4400 	1.4400 	1.5338 	1.6072 	1.8468 	1.8947  ];
DRO5500 = [3.6891	3.6586 	3.6586 	3.6586 	3.6586 	3.4902 	3.549626682	3.560208927	3.570791172];

 
RS500 = [1.4725	1.440 	1.5115 	1.5127 	1.5127 	1.5127 	1.5338 	1.6072 	1.6806 ];
RS5500 = [3.6891 3.6171 	3.5103 	3.4558 	3.4558 	3.4036 	3.3602 	3.5658 	3.7915 ];

x = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16]; 
x2 = [1 1.05 1.1 1.2 1.3 1.4 1.5 1.6 1.7] ;


figure()
hold on
dro5_5_aver = plot(x,DRO500, '-ko') ;
dro5_25_aver = plot(x, DRO5500,'-kv');
% 画出mu=0时的测试解，作为比较。并画出expost解，作为benchmark，看做最佳解。
l5 = plot(x,[1.42 1.42 1.42 1.42 1.42 1.42 1.42 1.42 1.42],'k-')
l55 = plot(x,[3.3168 3.3168  3.3168 3.3168 3.3168 3.3168 3.3168 3.3168 3.3168], 'k--')
 l155 = plot(0,1.4725,'s', 'MarkerFaceColor','k')
l255 =  plot(0,3.6891,'s', 'MarkerFaceColor','k')
  legend([dro5_5_aver dro5_25_aver l5 l55 l255], 'B_w=0.5km', 'B_w=5.5km', 'Ex post (B_w=0.5km)', 'Expost (B_w=5.5km)','Empirical');
 
xlabel('\mu')
ylabel('Average distance (km)')
grid on



figure()
hold on
 RS5_5_aver = plot(x2,RS500, '-ko' ) ;
 RS5_25_aver = plot(x2, RS5500,'-kv')
l5 = plot(x2,[1.42 1.42 1.42 1.42 1.42 1.42 1.42 1.42 1.42],'k-')
l55 = plot(x2,[3.3168 3.3168  3.3168 3.3168 3.3168 3.3168 3.3168 3.3168 3.3168], 'k--')
l55=  plot(1,1.4725,'s', 'MarkerFaceColor','k')
l225 =  plot(1,3.6891,'s', 'MarkerFaceColor','k')
% 画出tau=1时的测试解，作为比较。并画出expost解，作为benchmark，看做最佳解。
  
%   legend([RS5_5_aver RS5_25_aver l5 l55 l55], 'B_w=0.5km', 'B_w=2.5km', 'Ex post (B_w=0.5km)', 'Expost (B_w=2.5km)','Empirical');
 
 xlabel('\tau/Z_0')
ylabel('Average distance (km)')
grid on
xlim([1,1.6])
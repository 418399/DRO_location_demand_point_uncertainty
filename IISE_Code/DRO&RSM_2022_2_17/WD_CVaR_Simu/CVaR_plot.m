% Expose ;SAA;DRO;RS

MAX = [3.8503E+03	5.3916E+03	7.8005E+03	9.8000E+03	1.2021E+04	1.3000E+04
4.222E+03	5.999E+03	8.501E+03	1.120E+04	1.345E+04	1.500E+04
4.000E+03	5.610E+03	8.145E+03	1.024E+04	1.250E+04	1.400E+04
4.000E+03	5.600E+03	8.000E+03	9.991E+03	1.225E+04	1.355E+04]/1000 ;

MEAN = [1.437E+03	1.837E+03	2.358E+03	2.652E+03	2.936E+03	3.101E+03
1.505E+03	1.952E+03	2.529E+03	3.019E+03	3.522E+03	3.804E+03
1.505E+03	1.912E+03	2.434E+03	2.792E+03	3.253E+03	3.465E+03
1.527E+03	1.936E+03	2.440E+03	2.724E+03	3.137E+03	3.311E+03]/1000  ;

CVAR = [2.621E+03	3.609E+03	4.805E+03	5.544E+03	6.236E+03	6.389E+03
2.677E+03	3.789E+03	5.279E+03	6.511E+03	7.516E+03	8.190E+03
2.784E+03	3.746E+03	4.954E+03	5.878E+03	6.886E+03	7.465E+03
2.785E+03	3.73E+03	4.86E+03	5.68E+03	6.61E+03	7.17E+03] /1000 ;
 x = [0.5 1.5 2.5 3.5 4.5 5.5];

figure ()
hold on
max1 = plot(x, MAX(1,:), 'ko-') ;
max2 = plot(x, MAX(2,:), 'kv-') ;
max3 = plot(x, MAX(3,:), 'ks--') ;
max4 = plot(x, MAX(4,:), 'k^-') ;
 legend([max1 max2 max3 max4], 'EX post', 'SAA', 'C-DRO (\mu = 0.08)', 'C-RS (\tau/Z_0 = 1.4)' );
xlabel('Bandwidth (km)')
ylabel('Maximum distance (km)')
grid on

figure ()
hold on
cvar1 = plot(x, CVAR(1,:), 'ko-') ;
cvar2 = plot(x, CVAR(2,:), 'kv-') ;
cvar3 = plot(x, CVAR(3,:), 'ks--') ;
cvar4 = plot(x, CVAR(4,:), 'k^-') ;
% legend([cvar1 cvar2 cvar3 cvar4], 'EX post', 'SAA', 'C-DRO (\mu = 0.08)', 'C-RS (\tau/Z_0 = 1.4)' );
xlabel('Bandwidth (km)')
ylabel('CVaR distance (km)')
grid on



figure ()
hold on
mean1 = plot(x, MEAN(1,:), 'ko-') ;
mean2 = plot(x, MEAN(2,:), 'kv-') ;
mean3 = plot(x, MEAN(3,:), 'ks--') ;
mean4 = plot(x, MEAN(4,:), 'k^-') ;
xlabel('Bandwidth (km)')
ylabel('Mean distance (km)')
grid on



clc;clear; 
tic
load His ;
%load DC_20 ;
His = Sourcedata1(1:10,:);  % Historical data point
%[DC] = Distance(DC ) ;
b = 5;  % The max number of opened DC
aa = 5 ;  % number of demand points in x-axes
%bb = 4 ;  % ...in y-axes
  aaa = 4 ;
  bbb = 3 ;
multiplier = 0.3 ;
[His] = Distance(His) ;
% minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
% maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
 minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
lengthY = maxlo(2) - minlo(2) ;
EmProba = 1/size(His,1) ;
Disc_X = [minlo(1):lengthX/aa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_Y = [minlo(2):lengthX/aa:maxlo(2)] ;

LC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_X,2)
    for j = 1:size(Disc_Y,2)
        LC = [LC ; Disc_X(i) , Disc_Y(j) ] ;
    end
end

%% DC
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end

%%
 a = [];  % Find the MAX theta
% �ҳ������еļ���
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % �о���������귶Χ
EmProba = 1/size(His,1) ;  %ÿ����ʷ���ݵ�ȡ���ĸ���
temp = 0;
for i = 1:size(area,1)
    for j = 1:size(His,1)
        temp = temp + norm(area(i,:) - His(j,:),1) ;
    end
    a = [a;temp*EmProba];
    temp = 0;
end
%maxtheta = max(a) ;

maxtheta = max(a)  ; 
%+++++++++++++++++++++++++++++++++++++++
theta = multiplier * maxtheta ;  %  �������theta��0.1����Ϊ�뾶

%% ��ģ
iter = 1;
Sbar = LC ; % ��ʼ��z���ڵ��ɳ��˵�S������ɢ���ʾ
stopcri = 0.0005 ; % �㷨ֹͣ׼��
stop = 0;
%if_all_feasible = [];
LB = [] ;
UB = [];
New_Sbar = [];
New_His = [];
%tempAn = [];
%save_extremepoints = {};
%% ѭ�����Ĺ���
while stop==0
    %% 1������ɳ�����
    Sbar = [Sbar;New_Sbar];
    [xi,an,beta,Obj] = optimize(theta,b,His,DC,Sbar,EmProba);
    savesolutions{iter,1} = xi ; %��¼ÿ�ε����Ľ��Ŀ�꺯��ֵ
    savesolutions{iter,2} = an ;
    savesolutions{iter,3} = beta ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj] ;
    %% 2���жϵõ��Ľ��Ƿ���������Լ������
    % ����õĿ���DC��xi==1�����վ���z��������ֲ�ͬ������
    area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;
    Center = DC(find(xi==1),:) ; % �ҵ����ŵ�DC������
    Vcells = Voronoi(Center,area)  ;  %  ��ȡÿ��Pi�ı߽�㣨���㣩
    close()
    for n = 1:size(His,1) % ��ÿ��n \in N������֤�Ƿ����
        
        temp_Ix = [];
        for ii = 1:size(Vcells,1)  % Ѱ�Ҳ�ͬPi�£�Լ���Ҳ�����ֵ
            temp_Vn = [];
            xy = [Vcells{ii,1}  Vcells{ii,2}] ;  % ����S�Ķ�������
            xy = xy(1:end-1,:) ;  %ȥ����β��ͬ�ĵ�
            HisPoint = His(n,:);
            [extremepoints] = extreme(xy,HisPoint) ;  %  ����õ�z���������Ӧ�� w
            extremepoints = extremepoints(:,1:end-1) ;
            save_extremepoints{n,ii} = extremepoints;
            for jj = 1:size(extremepoints,1)
                bb = norm(Center(ii,:) - extremepoints(jj,1:2), 2 ) - beta * extremepoints(jj,3) ; % �������м��㣬����ai ��ÿ������֮��ľ���
                temp_Vn = [temp_Vn  bb];  %��¼���룬����ҵ����ֵ
            end
            [max_Vn Vn_z] = max(temp_Vn) ;  %�ҵ����м����У�1��2��������ֵ���Լ��ĸ�zʹ�����
            temp_Ix = [temp_Ix; max_Vn,extremepoints(Vn_z,1:2)] ; %��ÿ�� Pi�е����ֵ��¼����
        end   %end of ii
        
        [max_Ix Ix_i]= max(temp_Ix(:,1)) ; %�ҵ�����i�е����ֵ�Լ��ĸ�i����an_bar �Ƚ�
        % �ҵ� ʹ��Լ���Ҷ����� i��z��w��������룬����Ҷ����ֵ
        An(n) =  max_Ix ; %����Ҷ����ֵ
        %   save_An_z{n} = [An(n) temp_Ix(Ix_i,2:end)]   ;
        Zn(n,:) = temp_Ix(Ix_i,2:end) ;
        n
    end %  end of n \in N
    
    %  if_feasibleNum = an>=An'   %
    if_feasibleNum = roundn(an-An',-5) >=0 ;%
    tempAn = max(An,0) ;
    UB = [UB;  beta * theta + EmProba*sum(tempAn)  ] ; %�����������n������Լ������AnΥ��������õ��Ͻ�
    if sum(if_feasibleNum) < size(His,1)
        inde = find(if_feasibleNum==0) ;
        New_Sbar = [New_Sbar;  Zn(inde,:)] ;
    end
    
    
    if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)   || if_stop <= stopcri % ��������Լ������ﵽһ�����ȣ�ֹͣ�㷨
        stop = 1;
    end
    
    iter = iter + 1
    
    %     if iter>=20
    %         stop = 1 ;
    %     end
    
end %  end of while

x =  savesolutions{end,1}
Obj
toc


figure()
open = find(xi ==1)  ;
plot( DC(open,1),DC(open,2) , 'rs','MarkerFaceColor','k')
hold on
plot(LC(:,1), LC(:,2), 'k.')
axis equal
hold on
plot(His(:,1), His(:,2), 'bo')
hold on
plot (DC(:,1), DC(:,2), 's')






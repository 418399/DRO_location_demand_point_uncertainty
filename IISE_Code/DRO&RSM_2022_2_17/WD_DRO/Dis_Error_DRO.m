%%%%%%%%%%%
%% ������Ͻ��ѡȡ UB ��ѡȡ

clc;clear;
tic
load His ;
 His = Sourcedata1(1:10,:);  % Historical data point
%  [DC] = Distance(DC) ;
b = 5;  % The max number of opened DC
aa = 20;  % number of demand points in x-axes
%bb = 29   ;% ...in y-axes
  aaa = 4 ;
  bbb = 3 ;
 multiplier = 0.3;  %The factor to control the parameter \theta
[His] = Distance(His) ;
 minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
lengthY = maxlo(2) - minlo(2) ;

Disc_X = [minlo(1):aa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_Y = [minlo(2):aa:maxlo(2)] ;  size(Disc_X,2)*size(Disc_Y,2)
L_delta = aa ;
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

 
EmProba = 1/size(His,1) ;
%% get the max theta allover the region
a = [];  % Find the MAX theta
for j = 1:size(LC,1)
    for n = 1:size(His,1)
        a(j,n) = norm(His(n,:) - LC(j,:),1  ) ;
    end
end
a = sum(a,2) ;
maxtheta = max(a) * EmProba ; %  The max Wasserstein distance in the region
%+++++++++++++++++++++++++++++++++++++++++++

%theta = multiplier * maxtheta  
%% %
  theta = 4.681126646568384e+03;
 
%+++++++++++++++++++++++++++++++++++++++++++
No_initialLC = round(0.005* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;
%% ��ģ
iter = 1;
Sbar = LC(sort(random_Num'),:) ; % ��ʼ��z���ڵ��ɳ��˵�S������ɢ���ʾ
stopcri = 0.0005 ; % �㷨ֹͣ׼��
stop = 0;
LB = [] ;
UB = [];
New_Sbar = [];

for i=1:size(DC,1)  % DC������LC��2����
    for j = 1:size(LC,1)
        save_2_norm(i,j) = norm( DC(i,:) - LC(j,:) ,2) ;
    end
end
 
%% ѭ�����Ĺ���
while stop==0
    %% 1������ɳ�����
    Sbar = [Sbar;New_Sbar];
    Sbar = unique(Sbar,'rows')  ; %ȥ���ظ�����
    [xi,an,beta,Obj  ] = optimize(theta,b,His,DC,Sbar,EmProba);
    savesolutions{iter,1} = xi ; %��¼ÿ�ε����Ľ��Ŀ�꺯��ֵ
    savesolutions{iter,2} = an ;
    savesolutions{iter,3} = beta ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj] ;
    %% 2���жϵõ��Ľ��Ƿ���������Լ������
    
    for n =1:size(His,1)
        HisPoint = His(n,:);
        temp_save_2_norm = save_2_norm( find(xi==1),:)  ;
        [min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % ����ÿ��z���ҵ����������DC������
        
        for j = 1:size(LC,1)  % ���ݵõ���beta��kesi�õ� ���� N �µ�1����ֵ
            save_1_norm(n,j) = beta * norm( HisPoint - LC(j,:) ,1) ;
        end
        
        for j = 1:size(LC,1)  %
            temp(j) = min_value(j) - save_1_norm(n,j) ;
        end
        [max_An Lo_Z]= max(temp)  ; % ��ÿ��n�����еõ� ���ֵ An
        An(n) = max_An ;
        Zn(n,:) = LC(Lo_Z,:) ;
       
    end  % end of n
    

    if_feasibleNum = roundn(an-An',-5) >=0 ;%
 %   if_feasibleNum = an>=An'   %
    if_feasibleNum'
    tempAn = max(An,0) ;
    UB = [UB;  beta * theta + EmProba*sum(tempAn)  ] ; %�����������n������Լ������AnΥ��������õ��Ͻ�
    if sum(if_feasibleNum) < size(His,1)
        inde = find(if_feasibleNum==0) ;
        New_Sbar = [New_Sbar;  Zn(inde,:)] ;
    end
    
    
    if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)    %|| if_stop <= stopcri % ��������Լ����ֹͣ�㷨
        stop = 1;
    end
    
    iter = iter + 1
 
end
x =  savesolutions{end,1}
Obj
L_delta
 toc
 


plot(LC(:,1), LC(:,2), 'k.')
axis equal
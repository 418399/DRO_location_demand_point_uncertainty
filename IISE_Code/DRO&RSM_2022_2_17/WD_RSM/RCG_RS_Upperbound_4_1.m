clc;clear;
tic
load His ;
%  load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
% [DC] = Distance(DC) ;
b = 10;  % The max number of opened DC
aa = 39;  % number of demand points in x-axes
bb = 39  ;% ...in y-axes
MultiFactor = 1.05 ;
aaa = 7 ;
bbb = 4 ;

[His] = Distance(His) ;
% minlo = [min([His(:,1) ; DC(:,1)]) min([His(:,2) ;DC(:,2)])] ;
% maxlo = [max([His(:,1) ;DC(:,1)]) max([His(:,2);DC(:,2)])];
minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
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

EmProba = 1/size(His,1) ;
% �õ� Z_0
[z_0 x_0] = getZ0(EmProba, DC, His, b) ;
% ѡ�� tau  +++++++++++++++++++++++++

tau = z_0 * MultiFactor ;
No_initialLC = round(0.005* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;

% ��ģ
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



while stop==0
    % 1������ɳ�����
    Sbar = [Sbar;New_Sbar];
    Sbar = unique(Sbar,'rows')  ; %ȥ���ظ�����
    [xi,yij, sn, k, Obj] = optimize_RSM(EmProba, DC, Sbar, His, b, tau)  ;
    savesolutions{iter,1} = xi ; %��¼ÿ�ε����Ľ��Ŀ�꺯��ֵ
    savesolutions{iter,2} = sn ;
    savesolutions{iter,3} = k ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj]
    % 2���жϵõ��Ľ��Ƿ���������Լ������
    for n =1:size(His,1)
        HisPoint = His(n,:);
        temp_save_2_norm = save_2_norm( find(xi==1),:)  ;
        [min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % ����ÿ��z���ҵ����������DC������
        
        for j = 1:size(LC,1)  % ���ݵõ���beta��kesi�õ� ���� N �µ�1����ֵ
            save_1_norm(j) = norm( HisPoint - LC(j,:) ,1) ;
        end
        
        for j = 1:size(LC,1)  %���� k>=() �� ()�е���ҵ�()����֪ yn xi �µ����ֵ
            %  temp(j) = min_value(j) - save_1_norm(n,j) ;
            temp(j) = (min_value(j) - sn(n) ) / save_1_norm(j) ;
        end
        [max_An Lo_Z]= max(temp)  ; % ��ÿ��n�����еõ� ���ֵ An�Լ���Ӧ��zj
        tempKn(n,1) = max_An ;
        Zn(n,:) = LC(Lo_Z,:) ;
        
    end   % end of n
    Kn = max(tempKn) ;  %%%ͨ��max{ x/y }���k�����ֵ���Դ˷�������Sn, ������Ƿ�����  <=tau
    % �ж��Ƿ�������� zj ���ɳڽⶼ����Լ������yn_bar�������ɳڽ⣬���� sum(yn_bar)*1/N  <= tau��
    % ���ǣ����ɳڽ�Ϊ���н⣬��Ϊ���Ž�
    
    %% �ж��Ƿ�����Լ���� sn>=.....
    for nn = 1:size(His,1)
        Selected_Xi = DC(find(xi ==1 ),:)  ;  % ����ѡ�е�zj �����ŵ� DC֮��ľ���
        for ii=1:size(Selected_Xi,1)  % DC������LC��2����
            for jj = 1:size(Zn,1)  %
                zj_distance(ii,jj) = norm( Selected_Xi(ii,:) - Zn(jj,:) ,2) ;
            end
        end
        [Sub_min_Distance Sub_index_zj  ] = min(zj_distance,[],1)  ;
        Jus_sn(nn) = Sub_min_Distance(nn)  - k * norm( His(nn,:) - Zn(nn,:) ,1) ; % ���κ�zj������An
    end
    if_feasibleNum = roundn(sn-Jus_sn',-5) >=0 ;% �ж��Ƿ�����Լ������Ϊ1��������Լ��
    New_Sbar = [New_Sbar;  Zn(find(if_feasibleNum==0),:)] ;
    
%% �Ͻ磺
for nnn = 1:size(His,1)
 tempRn(nnn) = ( Sub_min_Distance(nnn)-sn(nnn) ) / norm( His(nnn,:) - Zn(nnn,:) ,1) ;
end
worst_R = max(tempRn)  ;
 
  UB = [UB; min([UB; worst_R])]
    
    
    
    
    
    
    
    
%     if_feasibleNum = roundn(k-Kn,-5) >=0 ;% �ж��Ƿ�����Լ������Ϊ1��������Լ��
%     New_Sbar = [New_Sbar;  Zn(find(if_feasibleNum==0),:)] ;
%     
%     
%     %����1���� �ҵ��� zj ����ӽ� �Ӽ��� ѭ����⡣
%     %            2������k_bar, ���� k_bar �� zj_bar,  �õ� yn_new�� ��֤�Ƿ�С�� tau��
%     %            ���ǣ���Ϊ���н⣬�õ��Ͻ� k_bar_new
%     if sum(if_feasibleNum) < size(His,1)
%         %% ����ѡ�е�zj �����ŵ� DC֮��ľ���
%         Selected_Xi = DC(find(xi ==1 ),:)  ;
%         for i=1:size(Selected_Xi,1)  % DC������LC��2����
%             for j = 1:size(Zn,1)  %
%                 zj_distance(i,j) = norm( Selected_Xi(i,:) - Zn(j,:) ,2) ;
%             end
%         end
%         [Sub_min_Distance Sub_index_zj  ] = min(zj_distance,[],1)  ;   % ����n���õ���zj���ҵ����������DC������
%         not_sati = [];  %��¼������Լ����zj������
%         for j = 1:size(Zn,1)  %��֤ÿ��ʹ��ȡ��max�� bar z_j �������� sn ������µ�Kn, �Լ��ɳڽ⣨x_i������ɵĽ��Ƿ�����Լ����
%             for n = 1:size(His,1)
%                 New_An(n) =  Sub_min_Distance(j) - Kn(n) * norm( His(n,:) - Zn(j,:) ,1) ; % ���κ�zj������An
%             end
%             if EmProba * sum(New_An) <= tau
%                 not_sati = [not_sati 0];
%             else
%                 not_sati = [not_sati 1];
%             end
%         end
%         if sum(not_sati) == size(Zn,1)  % ������е�zj�������� <= tau, ���Ͻ粻����
%             Ub = tau * 10e6 ;
%         else
%             Ub  = min(   Kn(find(not_sati==0))  ) ;
%         end
%         
%         UB = [UB; Ub]
%     end
    
        if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)     || if_stop <= stopcri % ��������Լ����ֹͣ�㷨
        stop = 1;
    end
    if_feasibleNum'
    iter = iter + 1
    
end
x =  savesolutions{end,1}
Obj
toc
%   [Obj_Value1] = Func_Test_RSM0001(x')
%    [Obj_Value5] = Func_Test_RSM0005(x')




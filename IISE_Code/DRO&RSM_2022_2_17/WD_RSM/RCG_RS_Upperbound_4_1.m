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
% 得到 Z_0
[z_0 x_0] = getZ0(EmProba, DC, His, b) ;
% 选择 tau  +++++++++++++++++++++++++

tau = z_0 * MultiFactor ;
No_initialLC = round(0.005* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;

% 建模
iter = 1;
Sbar = LC(sort(random_Num'),:) ; % 初始化z所在的松弛了的S区域，离散点表示
stopcri = 0.0005 ; % 算法停止准则
stop = 0;
LB = [] ;
UB = [];
New_Sbar = [];

for i=1:size(DC,1)  % DC到所有LC的2范数
    for j = 1:size(LC,1)
        save_2_norm(i,j) = norm( DC(i,:) - LC(j,:) ,2) ;
    end
end



while stop==0
    % 1，求解松弛问题
    Sbar = [Sbar;New_Sbar];
    Sbar = unique(Sbar,'rows')  ; %去除重复的行
    [xi,yij, sn, k, Obj] = optimize_RSM(EmProba, DC, Sbar, His, b, tau)  ;
    savesolutions{iter,1} = xi ; %记录每次迭代的解和目标函数值
    savesolutions{iter,2} = sn ;
    savesolutions{iter,3} = k ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj]
    % 2，判断得到的解是否满足所有约束条件
    for n =1:size(His,1)
        HisPoint = His(n,:);
        temp_save_2_norm = save_2_norm( find(xi==1),:)  ;
        [min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % 对于每个z，找到与其最近的DC及距离
        
        for j = 1:size(LC,1)  % 根据得到的beta和kesi得到 任意 N 下的1范数值
            save_1_norm(j) = norm( HisPoint - LC(j,:) ,1) ;
        end
        
        for j = 1:size(LC,1)  %计算 k>=() 中 ()中的项，找到()在已知 yn xi 下的最大值
            %  temp(j) = min_value(j) - save_1_norm(n,j) ;
            temp(j) = (min_value(j) - sn(n) ) / save_1_norm(j) ;
        end
        [max_An Lo_Z]= max(temp)  ; % 在每次n迭代中得到 最大值 An以及对应的zj
        tempKn(n,1) = max_An ;
        Zn(n,:) = LC(Lo_Z,:) ;
        
    end   % end of n
    Kn = max(tempKn) ;  %%%通过max{ x/y }求得k的最大值，以此反过来求Sn, 并检查是否满足  <=tau
    % 判断是否对于所有 zj ，松弛解都满足约束。（yn_bar由于是松弛解，满足 sum(yn_bar)*1/N  <= tau）
    % 若是，则松弛解为可行解，即为最优解
    
    %% 判断是否满足约束， sn>=.....
    for nn = 1:size(His,1)
        Selected_Xi = DC(find(xi ==1 ),:)  ;  % 计算选中的zj 到开放的 DC之间的距离
        for ii=1:size(Selected_Xi,1)  % DC到所有LC的2范数
            for jj = 1:size(Zn,1)  %
                zj_distance(ii,jj) = norm( Selected_Xi(ii,:) - Zn(jj,:) ,2) ;
            end
        end
        [Sub_min_Distance Sub_index_zj  ] = min(zj_distance,[],1)  ;
        Jus_sn(nn) = Sub_min_Distance(nn)  - k * norm( His(nn,:) - Zn(nn,:) ,1) ; % 对任何zj，计算An
    end
    if_feasibleNum = roundn(sn-Jus_sn',-5) >=0 ;% 判断是否都满足约束，若为1，则满足约束
    New_Sbar = [New_Sbar;  Zn(find(if_feasibleNum==0),:)] ;
    
%% 上界：
for nnn = 1:size(His,1)
 tempRn(nnn) = ( Sub_min_Distance(nnn)-sn(nnn) ) / norm( His(nnn,:) - Zn(nnn,:) ,1) ;
end
worst_R = max(tempRn)  ;
 
  UB = [UB; min([UB; worst_R])]
    
    
    
    
    
    
    
    
%     if_feasibleNum = roundn(k-Kn,-5) >=0 ;% 判断是否都满足约束，若为1，则满足约束
%     New_Sbar = [New_Sbar;  Zn(find(if_feasibleNum==0),:)] ;
%     
%     
%     %若否，1、则将 找到的 zj 都添加进 子集中 循环求解。
%     %            2、更新k_bar, 根据 k_bar 和 zj_bar,  得到 yn_new， 验证是否小于 tau，
%     %            若是，则为可行解，得到上界 k_bar_new
%     if sum(if_feasibleNum) < size(His,1)
%         %% 计算选中的zj 到开放的 DC之间的距离
%         Selected_Xi = DC(find(xi ==1 ),:)  ;
%         for i=1:size(Selected_Xi,1)  % DC到所有LC的2范数
%             for j = 1:size(Zn,1)  %
%                 zj_distance(i,j) = norm( Selected_Xi(i,:) - Zn(j,:) ,2) ;
%             end
%         end
%         [Sub_min_Distance Sub_index_zj  ] = min(zj_distance,[],1)  ;   % 对于n个得到的zj，找到与其最近的DC及距离
%         not_sati = [];  %记录不满足约束的zj的数量
%         for j = 1:size(Zn,1)  %验证每个使得取到max的 bar z_j 所产生的 sn ，结合新的Kn, 以及松弛解（x_i）所组成的解是否满足约束。
%             for n = 1:size(His,1)
%                 New_An(n) =  Sub_min_Distance(j) - Kn(n) * norm( His(n,:) - Zn(j,:) ,1) ; % 对任何zj，计算An
%             end
%             if EmProba * sum(New_An) <= tau
%                 not_sati = [not_sati 0];
%             else
%                 not_sati = [not_sati 1];
%             end
%         end
%         if sum(not_sati) == size(Zn,1)  % 如果所有的zj都不满足 <= tau, 则上界不更新
%             Ub = tau * 10e6 ;
%         else
%             Ub  = min(   Kn(find(not_sati==0))  ) ;
%         end
%         
%         UB = [UB; Ub]
%     end
    
        if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)     || if_stop <= stopcri % 满足所有约束，停止算法
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




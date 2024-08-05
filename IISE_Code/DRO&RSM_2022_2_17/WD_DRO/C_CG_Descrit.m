clc;clear;
tic
load His ;
 %  load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
%  [DC] = Distance(DC) ;
b = 10;  % The max number of opened DC
aa = 29;  % number of demand points in x-axes
bb = 29   ;% ...in y-axes
  aaa = 7 ;
  bbb = 4 ;
%   m = [0.14 0.18 0.2 0.24 0.28] ;
%   save_solu_and_time = [];
%   for mmm =1:size(m,2)
multiplier = 0.18;  %The factor to control the parameter \theta

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

theta = multiplier * maxtheta ;

 
%+++++++++++++++++++++++++++++++++++++++++++
No_initialLC = round(1* size(LC,1) ) ;
random_Num = randperm(size(LC,1),No_initialLC) ;
%% 建模
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
 
%% 循环求解的过程
while stop==0
    %% 1，求解松弛问题
    Sbar = [Sbar;New_Sbar];
    Sbar = unique(Sbar,'rows')  ; %去除重复的行
    [xi,an,beta,Obj  ] = optimize(theta,b,His,DC,Sbar,EmProba);
    savesolutions{iter,1} = xi ; %记录每次迭代的解和目标函数值
    savesolutions{iter,2} = an ;
    savesolutions{iter,3} = beta ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj] ;
    %% 2，判断得到的解是否满足所有约束条件
    
    for n =1:size(His,1)
        HisPoint = His(n,:);
        temp_save_2_norm = save_2_norm( find(xi==1),:)  ;
        [min_value min_loca] = min(temp_save_2_norm,[],1)  ;   % 对于每个z，找到与其最近的DC及距离
        
        for j = 1:size(LC,1)  % 根据得到的beta和kesi得到 任意 N 下的1范数值
            save_1_norm(n,j) = beta * norm( HisPoint - LC(j,:) ,1) ;
        end
        
        for j = 1:size(LC,1)  %
            temp(j) = min_value(j) - save_1_norm(n,j) ;
        end
        [max_An Lo_Z]= max(temp)  ; % 在每次n迭代中得到 最大值 An
        An(n) = max_An ;
        Zn(n,:) = LC(Lo_Z,:) ;
       
    end  % end of n
    

    if_feasibleNum = roundn(an-An',-5) >=0 ;%
 %   if_feasibleNum = an>=An'   %
    if_feasibleNum'
    tempAn = max(An,0) ;
    UB = [UB;  beta * theta + EmProba*sum(tempAn)  ] ; %如果不是所有n都满足约束，则An违反，带入得到上界
    if sum(if_feasibleNum) < size(His,1)
        inde = find(if_feasibleNum==0) ;
        New_Sbar = [New_Sbar;  Zn(inde,:)] ;
    end
    
    
    if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)    %|| if_stop <= stopcri % 满足所有约束，停止算法
        stop = 1;
    end
    
    iter = iter + 1
 
end
x =  savesolutions{end,1}
Obj
time = toc
 
% save_solu_and_time = [save_solu_and_time [x;Obj;time]]
%   end

%   [Obj_Value1] = Func_Test_RSM0001(x')
%    [Obj_Value5] = Func_Test_RSM0005(x')
% figure()
% open = find(xi ==1)  ;
% plot( DC(open,1),DC(open,2) , 'rs','MarkerFaceColor','k')
% hold on
% plot(LC(:,1), LC(:,2), 'k.')
% axis equal
% hold on
% plot(His(:,1), His(:,2), 'bo')
% hold on
% plot (DC(:,1), DC(:,2), 's')



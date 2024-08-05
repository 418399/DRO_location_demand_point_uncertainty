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
% 找出区域中的极点
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % 研究区域的坐标范围
EmProba = 1/size(His,1) ;  %每个历史数据点取到的概率
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
theta = multiplier * maxtheta ;  %  设置最大theta的0.1倍作为半径

%% 建模
iter = 1;
Sbar = LC ; % 初始化z所在的松弛了的S区域，离散点表示
stopcri = 0.0005 ; % 算法停止准则
stop = 0;
%if_all_feasible = [];
LB = [] ;
UB = [];
New_Sbar = [];
New_His = [];
%tempAn = [];
%save_extremepoints = {};
%% 循环求解的过程
while stop==0
    %% 1，求解松弛问题
    Sbar = [Sbar;New_Sbar];
    [xi,an,beta,Obj] = optimize(theta,b,His,DC,Sbar,EmProba);
    savesolutions{iter,1} = xi ; %记录每次迭代的解和目标函数值
    savesolutions{iter,2} = an ;
    savesolutions{iter,3} = beta ;
    savesolutions{iter,4} = Obj;
    savesolutions{iter,5} = New_Sbar ;
    New_Sbar = [];
    LB = [LB; Obj] ;
    %% 2，判断得到的解是否满足所有约束条件
    % 将求得的开放DC，xi==1，按照距离z最近，划分不同的区域
    area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;
    Center = DC(find(xi==1),:) ; % 找到开放的DC的坐标
    Vcells = Voronoi(Center,area)  ;  %  提取每个Pi的边界点（极点）
    close()
    for n = 1:size(His,1) % 对每个n \in N，都验证是否可行
        
        temp_Ix = [];
        for ii = 1:size(Vcells,1)  % 寻找不同Pi下，约束右侧的最大值
            temp_Vn = [];
            xy = [Vcells{ii,1}  Vcells{ii,2}] ;  % 区域S的顶点坐标
            xy = xy(1:end-1,:) ;  %去除首尾相同的点
            HisPoint = His(n,:);
            [extremepoints] = extreme(xy,HisPoint) ;  %  输出得到z的坐标和相应的 w
            extremepoints = extremepoints(:,1:end-1) ;
            save_extremepoints{n,ii} = extremepoints;
            for jj = 1:size(extremepoints,1)
                bb = norm(Center(ii,:) - extremepoints(jj,1:2), 2 ) - beta * extremepoints(jj,3) ; % 遍历所有极点，计算ai 与每个极点之间的距离
                temp_Vn = [temp_Vn  bb];  %记录距离，最后找到最大值
            end
            [max_Vn Vn_z] = max(temp_Vn) ;  %找到所有极点中，1、2范数最大的值，以及哪个z使其最大
            temp_Ix = [temp_Ix; max_Vn,extremepoints(Vn_z,1:2)] ; %将每个 Pi中的最大值记录下来
        end   %end of ii
        
        [max_Ix Ix_i]= max(temp_Ix(:,1)) ; %找到所有i中的最大值以及哪个i，与an_bar 比较
        % 找到 使得约束右端最大的 i，z，w，将其带入，获得右端最大值
        An(n) =  max_Ix ; %获得右端最大值
        %   save_An_z{n} = [An(n) temp_Ix(Ix_i,2:end)]   ;
        Zn(n,:) = temp_Ix(Ix_i,2:end) ;
        n
    end %  end of n \in N
    
    %  if_feasibleNum = an>=An'   %
    if_feasibleNum = roundn(an-An',-5) >=0 ;%
    tempAn = max(An,0) ;
    UB = [UB;  beta * theta + EmProba*sum(tempAn)  ] ; %如果不是所有n都满足约束，则An违反，带入得到上界
    if sum(if_feasibleNum) < size(His,1)
        inde = find(if_feasibleNum==0) ;
        New_Sbar = [New_Sbar;  Zn(inde,:)] ;
    end
    
    
    if_stop = (UB(end) - LB(end) ) / UB(end)
    if sum(if_feasibleNum) == size(His,1)   || if_stop <= stopcri % 满足所有约束，或达到一定精度，停止算法
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






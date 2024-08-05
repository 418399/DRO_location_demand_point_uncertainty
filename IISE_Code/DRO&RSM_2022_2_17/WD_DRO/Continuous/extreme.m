function [extremepoints] = extreme(xy,HisPoint)
%  xy = [a1 a2] ;  % 区域S的顶点坐标
% Hispoint   % 历史数据点
mt = 0; % z的自由度初始化为0
%% region中顶点的所有k子集
pointsw = [ ] ;
pointsz = [ ];
while mt <=2
    Pbar = nchoosek(1:2, mt ) ;
    Qbar = nchoosek(1:size(xy,1),size(xy,1) - mt - 1) ;
    if mt ==0
        tempw = [];
        for i = 1:size(xy,1)
            tempw = [tempw; norm(xy(i,:) - HisPoint ,1)] ; % w每次循环记录
        end
        % points{mt+1} = {tempw xy} ; % 将w，z存储到cell中
        pointsw = [pointsw; tempw] ;
        pointsz = [pointsz; xy];
    end
    
    %%
    if mt>0
        F = [];
        z = sdpvar(1,2,'full');
        lamda = sdpvar(size(xy,1),1,'full') ;
        for i = 1:size(Pbar,1)
            
            for k = 1:size(Pbar(i,:),2)
                F = [F,  z(k) == HisPoint(1,k) ];;  % Line 13 in the algorithm
            end
            
            for j = 1:size(Qbar,1)
                F = [F, z(1) == sum(lamda.*xy(:,1)  )  ];
                F = [F, z(2) ==  sum(lamda.*xy(:,2)  ) ];
                F = [F, sum(lamda) ==1 ];
                for kk = 1:size(Qbar,2)
                    F = [F, lamda(Qbar(j,kk)) == 0];  % 在每个Qbar子集下，把lamda==0写入约束
                end
                
                Obj = 0;
                ops = sdpsettings('solver','cplex','verbose',2 );
                [model,recoverymodel,diagnostic,internalmodel] = export(F,Obj,ops)  ;
                cons_ele = full(model.Aeq)  ;
                b_righthand = full(model.beq) ;
                RES=rref(cons_ele) ;  %化简行列式，看是否最优一行为0，即是否是线性无关的
                if sum(RES(end,:)) ==0
                    continue ;
                end
                answer = linsolve(cons_ele,b_righthand) ;  %线性无关则求解B（line15中的内容），得到z和lamda
                answer = roundn(answer,-3) ;  %设置解显示的精确度
                if isempty(find(answer(3:end)<0))==1  % 如果lamda全部为非负，则添加新极点
                    tw1 = norm( answer(1:2) - HisPoint' ,1 );  %添加w，（line18中的内容）
                    % points{mt+1} = {tw1  answer(1:2)'} ; % 将w，z存储到cell中
                    pointsw = [pointsw; tw1] ;
                    pointsz = [pointsz; answer(1:2)'];
                end
            end
        end   % end i
    end   % end mt>0
    mt = mt+1 ;
    
    extremepoints = [pointsz pointsw] ;  % 将z与w合并为一个数组
    
end   % end while

%% 去除非极点的点
F = [] ;
extremepoints(:,end+1) = 1:size(extremepoints,1)'  ;
for i = 1:size(extremepoints,1)
    tep = extremepoints;
   
    tep( find(tep(:,end)==i),:) = [] ;
 %   tep(i,:) = [];
    alpha = sdpvar(size(tep,1),1,'full');
    temp = 0;
    for j = 1:size(tep,1)
        temp = temp + alpha(j) * tep(j,1:end-1) ;
    end
    F = [extremepoints(find( extremepoints(:,end)==i ),1:end-1) == temp ] ;
    
    F = [F,alpha>=0];
    Obj = 0;
    ops = sdpsettings('solver','cplex','verbose',0 );
    solmp=solvesdp(F,Obj,ops );
    % solmp.info;
    alpha = value(alpha) ;
    if ~isempty(strfind(solmp.info,'Successfully'))
        extremepoints(find( extremepoints(:,end)==i ),:) = [];
    end
end
 
end


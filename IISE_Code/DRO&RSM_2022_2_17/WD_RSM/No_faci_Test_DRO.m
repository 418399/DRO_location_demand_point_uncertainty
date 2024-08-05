 clc;clear;
 %    load CDRO_No_facility1;
   load CRS_No_facility;;
   
   
 load Gene_data_His20_Band10000;
 for ii = 1:size(savex,2)
     xx = savex{ii}(1:end-1)' ;
 
aaa = 7 ;
bbb = 4 ;
 
load His ;
%   load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
[His] = Distance(His) ;
minlo = [min([His(:,1) ]) min([His(:,2) ])] ;
maxlo = [max([His(:,1) ]) max([His(:,2)])];
 
lengthX = maxlo(1) - minlo(1)  ;  % get the length of gap in x and y axes
lengthY = maxlo(2) - minlo(2) ;

 
Disc_XDC = [minlo(1):lengthX/aaa:maxlo(1)]  ;  % get the elements of the x and y  axes.
Disc_YDC = [minlo(2):lengthY/bbb:maxlo(2)] ;
DC = [];   %  Get the square lattice  locations
for i = 1:size(Disc_XDC,2)
    for j = 1:size(Disc_YDC,2)
        DC = [DC ; Disc_XDC(i) , Disc_YDC(j) ] ;
    end
end

%% 加工得到的测试点
Test_loca = [];
for i = 1:size(Gene_data,1)
    %   ttemp = Distance(Gene_data{i})  ;
    Test_loca = [Test_loca; Gene_data{i} ];
end
 
    %%  验证选址决策
    Selec_DC = DC(find(xx==1),:) ;
    temp = 0;
    for i = 1:size(Selec_DC,1)
        for j = 1:size(Test_loca,1)
            DistanceMat(i,j) = norm(Selec_DC(i,:) - Test_loca(j,:)  ,2) ;
        end
    end
    [Min_Dist_DC_Test Sub_index_x  ] = min(DistanceMat,[],1)  ;   % 对于n个得到的zj，找到与其最近的DC及距离
    temp = sum(Min_Dist_DC_Test)  ;
    Obj_Value(ii,1) = temp/100000  ;
 
 end

drox = [ 4 6 8 10 12 14 16] ;
plot(drox, Obj_Value, '-')







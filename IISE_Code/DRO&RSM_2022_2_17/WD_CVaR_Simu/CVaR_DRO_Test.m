function [averagedis] = CVaR_DRO_Test(x,DC,His) %将得到的x的解带入四川测试案例中，得到设施点到测试需求点的平均距离。
 
 
Selec_DC = DC(find(x==1),:) ;
temp = 0;
for i = 1:size(Selec_DC,1)
    for j = 1:size(His,1)
      %  DistanceMat(i,j) = distance(Selec_DC(i,1),Selec_DC(i,2),Test_loca(j,1),Test_loca(j,2))/180*pi*6371  ;  % km
         DistanceMat(i,j) = norm( DC(i,:) -Test_poi(j,:) ,2)  ;
    end
end


[Min_Dist_DC_Test Sub_index_x  ] = min(DistanceMat,[],1)  ;   % 对于n个得到的zj，找到与其最近的DC及距离
averagedis = sum(Min_Dist_DC_Test)/size(Min_Dist_DC_Test,2)  ;



end


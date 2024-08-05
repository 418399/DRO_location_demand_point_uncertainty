function [averagedis] = sichuantest(x,Sichuan_DC) %���õ���x�Ľ�����Ĵ����԰����У��õ���ʩ�㵽����������ƽ�����롣
Test_loca = [101.85	32.24
104.82	29.59
104.79	29.55
104.77	28.43
104.89	28.37
104.86	28.2
105	32.27
103	30.3
102.9	30.1
99.4	31
105.3	32.6
100.9	31.3
105.7	30.3
105.6	32.7
105.3	32.5
103.4	31.2
104.5	31.7
103.5	31.3
103.9	31.5
104.1	31.3];

Selec_DC = Sichuan_DC(find(x==1),:) ;
    temp = 0;
    for i = 1:size(Selec_DC,1)
        for j = 1:size(Test_loca,1)
            DistanceMat(i,j) = distance(Selec_DC(i,1),Selec_DC(i,2),Test_loca(j,1),Test_loca(j,2))/180*pi*6371  ;  % km
        end
    end


    [Min_Dist_DC_Test Sub_index_x  ] = min(DistanceMat,[],1)  ;   % ����n���õ���zj���ҵ����������DC������
    averagedis = sum(Min_Dist_DC_Test)/size(Min_Dist_DC_Test,2)  ;



end


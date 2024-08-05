clc;clear;
tic
load His ;
%  load DC_20 ;
His = Sourcedata1(1:20,:);  % Historical data point
%  [DC] = Distance(DC) ;
b = 3;   % The max number of opened DC
aa = 9;  % number of demand points in x-axes
bb = 9  ;% ...in y-axes
MultiFactor = 1.1;
aaa = 3 ;
bbb = 2.8006 ;
K=3;
 
mina = 0;  % define the initial values of the Bi-Section algorithm
maxb = 5;

%[His] = Distance(His) ;
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
[z_0 x_0] = getZ0(EmProba, DC, His, b) ;
tau = z_0 * MultiFactor ;


for i=1:size(DC,1)  % DC������LC��2����
    for j = 1:size(LC,1)
        save_2_norm(i,j) = norm( DC(i,:) - LC(j,:) ,2) ;
    end
end


load idx1;
load Center1 ;
His = [His idx ] ;  % add the cluster number of each points in the 3th column
%% Judge the candidate points belongs to which sub-region
area = [minlo(1) minlo(2); minlo(1) maxlo(2);  maxlo(1) maxlo(2);  maxlo(1) minlo(2)] ;  % �о���������귶Χ
Vcells = Voronoi(Center ,area)  ;
hold on
plot(His(:,1),His(:,2),'v')
%
close()
temp = 0;
for k = 1:K
    [in,on] = inpolygon(LC(:,1), LC(:,2), Vcells{k,1},Vcells{k,2});
    in_poly = double(in)  ;
    temp = temp + in_poly*k ;  % convert to the NO. of each demand point
end
LC(:,3) = temp ; % add the NO. to the 3th column

%% Didstinguish the clusters in different sub-regions
for k = 1:K
    DiffertKK(k,1) = { His( find(His(:,3) == k), 1:2 )  }  ;  % His
    DiffertKK(k,2) = {  LC( find( LC(:,3) == k ), 1:2  ) }  ; % LC
end

iter = 0;
stop=0;
while stop==0
    [Statu_b,xi_b] =  BiSectionConstraint(maxb,b,His,LC,DC,tau,DiffertKK,K) ; % ����֤rȡ�����ֵ��b�ˡ��������У���˵��ȡֵ����
    if Statu_b == 0
        Statu_b
        break;
    end
    
    [Statu_a, xi_a] =  BiSectionConstraint(mina,b,His,LC,DC, tau,DiffertKK,K) ; %����֤��Сֵa�Ƿ���У��������У�������
    if Statu_a == 1 %���aʱ���У�˵������Сֵʱ�ҵ����Ž⣬������Ž�xi
        xi_a
        mina
        break;
    end
    
    if Statu_a == 0 % ��aʱ�����У�������a��ֵ��ȡ(a+b)/2�����Ƿ��н�
        [Statu_c, xi_c] =  BiSectionConstraint((maxb + mina)/2,b,His,LC,DC, tau,DiffertKK,K) ; %����֤��Сֵa�Ƿ���У��������У�������
        if Statu_c==1 % ��ʱ����м�ֵc�н⣬��˵�������½�c��ֵ
            maxb = (maxb + mina)/2
        else
            mina = (maxb + mina)/2
        end
    end
    if (maxb - mina)/maxb <=0.0001
        stop=1;
        if Statu_a==1
            xi_a
            mina
        elseif Statu_c==1
            xi_c
            (maxb + mina)/2
        elseif Statu_b==1
            xi_b
            maxb
        end
    end
    iter = iter + 1
end


toc

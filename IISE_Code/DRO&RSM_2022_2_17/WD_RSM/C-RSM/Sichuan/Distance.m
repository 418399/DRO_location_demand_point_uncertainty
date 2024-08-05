function [new_position] = Distance(Matrix)
 
M_PI=3.14159265358979323846;
L = 6381372 * M_PI * 2; %�����ܳ�  
W = L; % ƽ��չ����x������ܳ�  
H = L / 2; % y��Լ�����ܳ�һ��  
mill = 2.3; % ����ͶӰ�е�һ����������Χ��Լ������2.3֮�� 
n=size(Matrix,1);

new_position=[];
for i =1:n
    lon=Matrix(i,1);
    lat=Matrix(i,2);
    x = lon * M_PI / 180; % �����ȴӶ���ת��Ϊ����  
    y = lat * M_PI / 180; %��γ�ȴӶ���ת��Ϊ����  
    y1 = 1.25 * log(tan(0.25 * M_PI + 0.4 * y)); % ����ͶӰ��ת��  
    % ����תΪʵ�ʾ���  
    dikaerX = (W / 2) + (W / (2 * M_PI)) * x ; %�ѿ�������x
    dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;%�ѿ�������y
    new_position(i,1)=dikaerX;
    new_position(i,2)=dikaerY;
 
end
 




end


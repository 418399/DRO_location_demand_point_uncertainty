function [new_position] = Distance(Matrix)
 
M_PI=3.14159265358979323846;
L = 6381372 * M_PI * 2; %地球周长  
W = L; % 平面展开后，x轴等于周长  
H = L / 2; % y轴约等于周长一半  
mill = 2.3; % 米勒投影中的一个常数，范围大约在正负2.3之间 
n=size(Matrix,1);

new_position=[];
for i =1:n
    lon=Matrix(i,1);
    lat=Matrix(i,2);
    x = lon * M_PI / 180; % 将经度从度数转换为弧度  
    y = lat * M_PI / 180; %将纬度从度数转换为弧度  
    y1 = 1.25 * log(tan(0.25 * M_PI + 0.4 * y)); % 米勒投影的转换  
    % 弧度转为实际距离  
    dikaerX = (W / 2) + (W / (2 * M_PI)) * x ; %笛卡尔坐标x
    dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;%笛卡尔坐标y
    new_position(i,1)=dikaerX;
    new_position(i,2)=dikaerY;
 
end
 




end


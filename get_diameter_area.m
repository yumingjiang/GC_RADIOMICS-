function [diameter,myArae] =get_diameter_area( mask )
%function get_diameter_area:获取mask的直径和面积
%description：
%输入mask，值不为0的地方为感兴趣区域，求其直径和面积
%INPUTS:
%mask:处理的矩阵
%OUTPUTS:
%diameter:直径
%area：面积
 
[m,n]=size(mask);
myArae=0;
diaTemp=[];
k=1;
for i=1:m
    for j=1:n
        if(mask(i,j)~=0)
            myArae=myArae+1;   %计算面积
            diaTemp(k,1)=i;  %区域存储坐标
            diaTemp(k,2)=j;
            k=k+1;
        else
            continue;
        end
    end
end
 
n=size(diaTemp);
length=[];
k=1;
for i=1:n
    for j=i+1:n
        length(k)=sqrt((diaTemp(j,1)-diaTemp(i,1)).^2+(diaTemp(j,2)-diaTemp(i,2)).^2);%计算任意两个坐标间的距离
        k=k+1;
    end
end
 
diameter=max(length);%得到最大距离，即直径
end
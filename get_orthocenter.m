 
function [orthocenter] = get_orthocenter( mask,area )
%function get_orthocenter:获取mask重心
%description：
%获取mask区域重心
%INPUTS:
%mask:需要计算的的图像
%OUTPUTS:
%orthocenter:重心
 
[m,n]=size(mask);
xTemp=0;
yTemp=0;
 
for i=1:m
    for j=1:n
        if(mask(i,j)~=0)
            xTemp=xTemp+i;
            yTemp=yTemp+j;
        else
            continue;
        end
    end
end
orthocenter.x=ceil(xTemp/area);
orthocenter.y=ceil(yTemp/area);
 
end
 

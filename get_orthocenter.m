 
function [orthocenter] = get_orthocenter( mask,area )
%function get_orthocenter:��ȡmask����
%description��
%��ȡmask��������
%INPUTS:
%mask:��Ҫ����ĵ�ͼ��
%OUTPUTS:
%orthocenter:����
 
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
 

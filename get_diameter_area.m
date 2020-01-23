function [diameter,myArae] =get_diameter_area( mask )
%function get_diameter_area:��ȡmask��ֱ�������
%description��
%����mask��ֵ��Ϊ0�ĵط�Ϊ����Ȥ��������ֱ�������
%INPUTS:
%mask:����ľ���
%OUTPUTS:
%diameter:ֱ��
%area�����
 
[m,n]=size(mask);
myArae=0;
diaTemp=[];
k=1;
for i=1:m
    for j=1:n
        if(mask(i,j)~=0)
            myArae=myArae+1;   %�������
            diaTemp(k,1)=i;  %����洢����
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
        length(k)=sqrt((diaTemp(j,1)-diaTemp(i,1)).^2+(diaTemp(j,2)-diaTemp(i,2)).^2);%�����������������ľ���
        k=k+1;
    end
end
 
diameter=max(length);%�õ������룬��ֱ��
end
function [shape] = get_shape(perimeter ,area )
%function get_shape:��ȡ��״����
%description��
%��ȡmask��״����
%INPUTS:
%perimeter ,area�����mask��Ҫ����ĵ�ͼ����ܳ������
%OUTPUTS:
%shape:��״����
 
shape=(perimeter.^2)/(4*pi*area);
end
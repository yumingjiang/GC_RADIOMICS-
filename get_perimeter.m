function [perimeter] = get_perimeter(mask)
%function get_perimeter :��ȡ�ܳ�
%purpose��
%��ȡmask�е��ܳ�
%INPUTS:
%mask��Ȧ������
%OUTPUTS:
%perimeter���ܳ�
 
m_edge=edge(mask);
edgeSpot=find(m_edge==1); %��ȡ�߽���������
perimeter=size(edgeSpot,1);
end
 


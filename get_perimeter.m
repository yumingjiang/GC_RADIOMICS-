function [perimeter] = get_perimeter(mask)
%function get_perimeter :获取周长
%purpose：
%获取mask中的周长
%INPUTS:
%mask：圈出区域
%OUTPUTS:
%perimeter：周长
 
m_edge=edge(mask);
edgeSpot=find(m_edge==1); %获取边界坐标数组
perimeter=size(edgeSpot,1);
end
 


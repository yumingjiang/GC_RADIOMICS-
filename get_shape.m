function [shape] = get_shape(perimeter ,area )
%function get_shape:获取形状参数
%description：
%获取mask形状参数
%INPUTS:
%perimeter ,area区域的mask需要计算的的图像的周长，面积
%OUTPUTS:
%shape:形状参数
 
shape=(perimeter.^2)/(4*pi*area);
end
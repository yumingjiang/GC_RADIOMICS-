function [CTmax,CTpeak,CTmean,aucCSH] = getCTmetrics(ROIonlyCT)
% -------------------------------------------------------------------------
% function [CTmax,CTpeak,CTmean,aucCSH] = getCTmetrics(ROIonlyCT)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes CTmax, CTpeak and CTmean, AUC-CSH and Percent 
% Inactive metrics from the region of interest (ROI) of an input CT volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonlyCT: 2D array representing the CT volume in CT format, with 
%               voxels outside the ROI set to NaNs. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - CTmax: Maximum SUV of the ROI.
% - CTpeak: Average of the voxel with maximum SUV within the ROI and its 
%            26 connected neighbours.
% - CTmean: Average CT value of the ROI.
% - aucCSH: Area under the curve of the cumulative CT-volume histogram
%           describing the percentage of total volume of the ROI above a 
%           percentage threshold of maximum CT.

ROIonlyCT = padarray(ROIonlyCT,[1 1 1],NaN);


[CTmax,indMax] = max(ROIonlyCT(:));


[indMaxX,indMaxY,indMaxZ] = ind2sub(size(ROIonlyCT),indMax);
connectivity = getneighbors(strel('arbitrary',conndef(3,'maximal')));
nPeak = length(connectivity);
neighborsMax = zeros(1,nPeak);
for i=1:nPeak
    neighborsMax(i) = ROIonlyCT(connectivity(i,1)+indMaxX,connectivity(i,2)+indMaxY,connectivity(i,3)+indMaxZ);
end
CTpeak = mean(neighborsMax(~isnan(neighborsMax)));

CTmean=mean(ROIonlyCT(~isnan(ROIonlyCT)));
[aucCSH] = getAUCCSH(ROIonlyCT);

end

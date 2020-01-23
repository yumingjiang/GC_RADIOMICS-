function Feature_Extr( dcm_path,nii_path, fund_path )

filename = fullfile(fund_path,'data.csv');
fid = fopen(filename,'w');

% ***************读取初始化************************************************
%%% dcm文件
I=dicomread(dcm_path); %读取图像

info_value = dicominfo(dcm_path);
pixel_spacing = info_value.PixelSpacing;
value_r_1 = 2/(pixel_spacing(1,:));
value_r_1 = round(value_r_1);
value_r_2 = 1/(pixel_spacing(1,:));
value_r_2 = round(value_r_2);

flag= min(min(I));
if flag < 0
    I = double(I);
    I = I + 2048;
    I = uint16(I);
end 

% I=dicomread('C:\Users\95720\Desktop\医学图像\roi1 540例\8\IMG-0001-00001.dcm'); %读取图像
% figure(1)
% imagesc(I);%显示图像
temp = double(I);
v= max(max(temp));
ss= min(min(temp));
s = temp*256/(v);
dcm_img = uint8(s);
% figure(2)
imshow(dcm_img)
imwrite(dcm_img,fullfile(fund_path,'dcm_img.jpg'));
% figure(3)
% imshow(I)
% img =uint16(I);              %将灰度级映射到0~255
% low=min(min(img));
% high=max(max(img));
% maxgray=high-low;           %计算窗宽
% rate=256/maxgray;
% img=I*rate;
% img=img+abs(min(min(img)));    %加窗
% dcm_img=uint8(img);           %转化为8位的位图数据格式
% figure(2)
% imshow(dcm_img)
% imwrite(dcm_img,fullfile(fund_path,'dcm_img.jpg'));



%%% nii文件
info  = load_nii (nii_path);
% info  = load_nii ('C:\Users\95720\Desktop\ROI 图像图\1\1.nii');
nii_img_r = info.img;

qq = size(nii_img_r,3); 
if qq > 1
    for i = 1:qq
        if max(max(nii_img_r(:,:,i))) == 1
            nii_img = nii_img_r(:,:,i);
            % disp('done!')
        else
            % continue
        end
    end
else
    nii_img = nii_img_r;
end
nii_img = imrotate(nii_img, 270,'nearest');
nii_img(:,:)=flipud(nii_img(:,:));

%% erode
se = strel('disk',value_r_1);
xxxx = imerode(nii_img(:,:),se);
%%
se1 = strel('disk',value_r_2);
yyyy = imdilate(nii_img(:,:),se1);
%%
[r,c] = find(xxxx==1);
num = length(r);

if num > 30
    nii_img(:,:) = yyyy-xxxx;
else
    disp(fund_path)
end
%%
% figure(3);
% imshow(nii_img,[])
% imwrite(nii_img,fullfile(fund_path,'nii_img.jpg'));
% *************************************************************************
% [ROIonly,levels] = prepareVolume(I,nii_img,'PETscan',4,3.27,1,5,'Matrix','Lloyd',32);



% *************原图和mask相与**********************************************
roi_img_w = dcm_img.*uint8(nii_img);
% figure(4)
% imshow(roi_img_w);
imwrite(roi_img_w,fullfile(fund_path,'roi_img_w.jpg'));
% *************************************************************************

% *****************************原图上显示roi*******************************
show_img = dcm_img;
[rows , cols ] = size(show_img);
%%% ***********最小拟合矩形
up_line = 255;
down_line = 0;
left_line = 255;
right_line = 0;
for i = 1:rows  
    for j = 1:cols  
        if nii_img(i,j) == 1
            if j >= right_line
                right_line = j;
            end
            if j <= left_line
                left_line = j;
            end
            if i <= up_line
                up_line = i;
            end
            if i >= down_line
                down_line = i;
            end
            show_img(i,j) = 255;
        end

    end 

end 
% figure(5)
% imshow(show_img)
imwrite(show_img,fullfile(fund_path,'show_img.jpg'));
% *************************************************************************

% **************************求出真的roi************************************
rect=[left_line up_line right_line-left_line down_line-up_line];%多边形的最小外接矩形 
roi_true=imcrop(dcm_img,rect);%裁剪ROI外接矩形，此处rect=[xmin ymin width hight] 
% imshow(roi_true);
imwrite(roi_true,fullfile(fund_path,'roi_true.jpg'));
% *************************************************************************

Hist=HistogramFeatures(double(dcm_img),nii_img);
cellperson = struct2cell(Hist);
len = length(cellperson);
for i=1:len
    fprintf(fid, '%07.7f\n', cellperson{i,1});
end    
fprintf(fid, '\n');
fprintf(fid, '\n');
    
    
shape = struct;
L = bwlabel(nii_img,8);
STATS = regionprops(L,'Area');
shape.area = STATS.Area;
STATS = regionprops(L,'Orientation');
shape.orientation = STATS.Orientation;
STATS = regionprops(L,'Eccentricity');
shape.eccentricity = STATS.Eccentricity;
STATS  = regionprops(L,'EquivDiameter');
shape.equivdiameter = STATS.EquivDiameter;
STATS   = regionprops(L,'Solidity');
shape.solidity = STATS.Solidity;
STATS  = regionprops(L,'Extent');
shape.extent = STATS.Extent;
STATS  = regionprops(L,'EulerNumber');
shape.eulernumber = STATS.EulerNumber;
STATS = regionprops(L,'perimeter');
shape.perimeter = STATS.Perimeter;

cellperson = struct2cell(shape);
len = length(cellperson);
for i=1:len
    fprintf(fid, '%07.7f\n', cellperson{i,1});
end    
fprintf(fid, '\n');
fprintf(fid, '\n');

% *************************************************************************
sigma_m = [1.0, 1.5,2.0, 2.5];
fprintf(fid, 'ORI\n');
for count = 1:5
    dcm_img_do = I;
    if count ~= 1
        %%%% 滤波
        fprintf(fid, '%1.1f滤波\n',sigma_m(count-1));
        sigma = sigma_m(count-1);
        gausFilter = fspecial('gaussian', [5,5], sigma);
        dcm_img_do = imfilter(dcm_img, gausFilter, 'replicate');
        % figure(2)
        % imshow(dcm_img)
        %%%% 
    end

    [ROIonly,levels] = prepareVolume(dcm_img_do,nii_img,'PETscan',4,3.27,1,5,'Matrix','Lloyd',32);
    % [ROIonly,levels] = prepareVolume(dcm_img_do,nii_img);
    [GLCM] = getGLCM(ROIonly,levels); %调用getGLCM获得GCLM矩阵
    % [glcmTextures] = getGLCMtextures(GLCM);%调用getGCLMtextures函数获得GCLM纹理
    [glcmTextures] = computeGLCMRadiomics(GLCM);%调用getGCLMtextures函数获得GCLM纹理 
    cellperson = struct2cell(glcmTextures);
    len = length(cellperson);
    % fprintf(fid, 'GLCM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');

    
    [GLRLM] = getGLRLM(ROIonly,levels); 
    [glrlmTextures] = getGLRLMtextures(GLRLM);
    cellperson = struct2cell(glrlmTextures);
    len = length(cellperson);
    % fprintf(fid, 'GLRLM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
    
    [GLSZM] = getGLSZM(ROIonly,levels); 
    [glszmTextures] = getGLSZMtextures(GLSZM);
    cellperson = struct2cell(glszmTextures);
    len = length(cellperson);
    % fprintf(fid, 'GLSZM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
    [NGTDM,countValid] = getNGTDM(ROIonly,levels); 
    [ngtdmTextures] = getNGTDMtextures(NGTDM,countValid);
    cellperson = struct2cell(ngtdmTextures);
    len = length(cellperson);
    % fprintf(fid, 'NGTDM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end







fclose(fid);

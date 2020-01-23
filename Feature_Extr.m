function Untitled( dcm_path,nii_path, fund_path )

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




% [globalTextures] = getGlobalTextures(ROIonly,100);
% [aucCSH] = getAUCCSH(ROIonly);
% [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(ROIonly);
% perimeter = get_perimeter(nii_img);
% [diameter,myarea] =get_diameter_area( nii_img );
% orthocenter = get_orthocenter( nii_img,myarea );
% shape = get_shape(perimeter ,myarea );
% Hu_Moment = HuSquare(roi_true);


%%%[G,gabout] = gaborfilter1(ROIonly,2,4,16,pi/3); 
%%% figure,imshow(gabout,[]);



% % % =================================================
% % 进行小波分解
% % % =================================================
% % 
% % 对图像进行2层的小波分解
% % N = 2 ;
% % [c,s] = wavedec2(ROIonly,N,'db1');
% % 
% % % =================================================
% % 提取各层对应的低频和高频系数
% % % =================================================
% % 
% % 对于原始图像512×512 ，其小波分解第1层维度为256×256，第2层维度为128×128
% % 提取小波分解中第1层低频系数和高频系数
% % a_ca1 = appcoef2(c,s,'db1',1);
% % a_ch1 = detcoef2('h',c,s,1);
% % a_cv1 = detcoef2('v',c,s,1);
% % a_cd1 = detcoef2('d',c,s,1);
% % 显示第1层的各分量
% % figure(6);
% % subplot(4,4,[3,4,7,8]);imshow(a_ch1,[]);
% % subplot(4,4,[9,10,13,14]);imshow(a_cv1,[]);
% % subplot(4,4,[11,12,15,16]);imshow(a_cd1,[]);
% % 
% % 提取第2层的低频系数和高频系数
% % ca2 = appcoef2(c,s,'db1',2);
% % ch2 = detcoef2('h',c,s,2);
% % cv2 = detcoef2('v',c,s,2);
% % cd2 = detcoef2('d',c,s,2);
% % % 显示第2层的各分量
% % subplot(4,4,1);imshow(ca2,[]);
% % subplot(4,4,2);imshow(ch2,[]);
% % subplot(4,4,5);imshow(cv2,[]);
% % subplot(4,4,6);imshow(cd2,[]);
% % 
% % 
% % % =================================================
% % 分别对各频率成分进行重构
% % % =================================================
% % 
% % 使用第2层进行重构
% % recon_a1 = wrcoef2('a',c,s,'db1',2);
% % recon_h1 = wrcoef2('h',c,s,'db1',2);
% % recon_v1 = wrcoef2('v',c,s,'db1',2);
% % recon_d1 = wrcoef2('d',c,s,'db1',2);
% % 显示各重构的成分 维度都在512×512
% % recon_set = [recon_a1,recon_h1;recon_v1,recon_d1];  
% % figure(7);imshow(recon_set,[]);title('第2层小波系数的重构');
% % 
% % % =================================================
% % 重构出原始图像
% % % =================================================
% % 
% % recon_img = recon_a1+recon_h1+recon_v1+recon_d1;
% % recon_img = mat2gray(recon_img );
% % figure(8);imshow(recon_img );title('重构出的原始图像');
% % WT = wavedec2(ROIonly, 2, 'db4');
% % LL = WT.dec;


fclose(fid);
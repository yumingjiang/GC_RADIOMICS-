function Untitled( dcm_path,nii_path, csv_path )

addpath('C:\Users\95720\Desktop\NIfTI_20140122')
filename = csv_path;
fid = fopen(filename,'a+');

% ***************读取初始化************************************************
%%% dcm文件
I=dicomread(dcm_path); %读取图像
% I=dicomread('C:\Users\95720\Desktop\ROI 图像图\1\ser004img00050.dcm'); %读取图像
% imagesc(I);%显示图像
img =double(I);              %将灰度级映射到0~255
low=min(min(img));
high=max(max(img));
maxgray=high-low;           %计算窗宽
rate=256/maxgray;
img=I*rate;
img=img+abs(min(min(img)));    %加窗
dcm_img=uint8(img);           %转化为8位的位图数据格式
% imshow(dcm_img)




%%% nii文件
info  = load_nii (nii_path);
% info  = load_nii ('C:\Users\95720\Desktop\ROI 图像图\1\1.nii');
nii_img = info.img;
nii_img = imrotate(nii_img, 270,'nearest');
nii_img(:,:)=flipud(nii_img(:,:));
% figure(3);
% imshow(nii_img,[])
% *************************************************************************



% *************原图和mask相与**********************************************
roi_img_w = dcm_img.*uint8(nii_img);
% figure(4)
% imshow(roi_img_w);
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
% *************************************************************************

% **************************求出真的roi************************************
rect=[left_line up_line right_line-left_line down_line-up_line];%多边形的最小外接矩形 
roi_true=imcrop(dcm_img,rect);%裁剪ROI外接矩形，此处rect=[xmin ymin width hight] 
% imshow(roi_true);
% *************************************************************************


% ******************求方差*************************************************
var_img = double(roi_true);
val_value=var(var_img(:)); 
% *************************************************************************

% *********************标准差**********************************************
std_img = double(roi_true);
std_value = std2(std_img);
% *************************************************************************

% **************************最大最小像素值及位置****************************
max_min_img = roi_true;
max_value = max(max(max_min_img));
[max_row max_col] = find(max_value == max_min_img);
min_value = min(min(max_min_img));
[min_row min_col] = find(min_value == max_min_img);
% *************************************************************************





% *******************************GLCM特征：能量、熵、惯性矩、相关***********
% GLCM = graycomatrix(roi_true,'GrayLimits',[]);
% stats = graycoprops(GLCM,{'energy','contrast','correlation','homogeneity'});

% gray = roi_true
% offsets = [0 1;-1 1;-1 0;-1 -1];
% m = 3; % 3阶灰度级
% [GLCMS,SI] = graycomatrix(gray,'GrayLimits',[],'Of',offsets,'NumLevels',m);
% P = GLCMS;
% [kk,ll,mm] = size(P);
% % 对共生矩阵归一化
% %---------------------------------------------------------
% for n = 1:mm
%    P(:,:,n) = P(:,:,n)/sum(sum(P(:,:,n)));
% end
%  
% %-----------------------------------------------------------
% %对共生矩阵计算能量、熵、惯性矩、相关4个纹理参数
% %-----------------------------------------------------------
% H = zeros(1,mm);
% I = H;
% Ux = H;      Uy = H;
% deltaX= H;  deltaY = H;
% C =H;
% for n = 1:mm
%    E(n) = sum(sum(P(:,:,n).^2)); %能量
%    for i = 1:kk
%       for j = 1:ll
%           if P(i,j,n)~=0
%              H(n) = -P(i,j,n)*log(P(i,j,n))+H(n); %熵
%           end
%           I(n) = (i-j)^2*P(i,j,n)+I(n);  %惯性矩
% 
%           Ux(n) = i*P(i,j,n)+Ux(n); %相关性中μx
%           Uy(n) = j*P(i,j,n)+Uy(n); %相关性中μy
%       end
%    end
% end
% for n = 1:mm
%    for i = 1:kk
%        for j = 1:ll
%             deltaX(n) = (i-Ux(n))^2*P(i,j,n)+deltaX(n); %相关性中σx
%             deltaY(n) = (j-Uy(n))^2*P(i,j,n)+deltaY(n); %相关性中σy
%             C(n) = i*j*P(i,j,n)+C(n);
%        end
%    end
%    C(n) = (C(n)-Ux(n)*Uy(n))/deltaX(n)/deltaY(n); %相关性
% end
%  
% %--------------------------------------------------------------------------
% %求能量、熵、惯性矩、相关的均值和标准差作为最终8维纹理特征
% %--------------------------------------------------------------------------
% a1 = mean(E)
% b1 = sqrt(cov(E))
%    
% a2 = mean(H)
% b2 = sqrt(cov(H))
%    
% a3 = mean(I)
% b3 = sqrt(cov(I))
%  
% a4 = mean(C)
% b4 = sqrt(cov(C))
%   
% sprintf('0,45,90,135方向上的能量依次为： %f, %f, %f, %f',E(1),E(2),E(3),E(4))  % 输出数据;
% sprintf('0,45,90,135方向上的熵依次为： %f, %f, %f, %f',H(1),H(2),H(3),H(4))  % 输出数据;
% sprintf('0,45,90,135方向上的惯性矩依次为： %f, %f, %f, %f',I(1),I(2),I(3),I(4))  % 输出数据;
% sprintf('0,45,90,135方向上的相关性依次为： %f, %f, %f, %f',C(1),C(2),C(3),C(4))  % 输出数据;
% *************************************************************************
sigma_m = [1.0, 1.5,2.0, 2.5];
fprintf(fid, 'ORI\n');
for count = 1:5
    dcm_img_do = dcm_img;
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
    
    [GLCM] = getGLCM(ROIonly,levels); %调用getGLCM获得GCLM矩阵
    [glcmTextures] = getGLCMtextures(GLCM);%调用getGCLMtextures函数获得GCLM纹理
    cellperson = struct2cell(glcmTextures);
    len = length(cellperson);
    fprintf(fid, 'GLCM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');

    
    [GLRLM] = getGLRLM(ROIonly,levels); 
    [glrlmTextures] = getGLRLMtextures(GLRLM);
    cellperson = struct2cell(glrlmTextures);
    len = length(cellperson);
    fprintf(fid, 'GLRLM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
    
    [GLSZM] = getGLSZM(ROIonly,levels); 
    [glszmTextures] = getGLSZMtextures(GLSZM);
    cellperson = struct2cell(glszmTextures);
    len = length(cellperson);
    fprintf(fid, 'GLSZM\n');
    for i=1:len
        fprintf(fid, '%07.7f\n', cellperson{i,1});
    end    
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
    [NGTDM,countValid] = getNGTDM(ROIonly,levels); 
    [ngtdmTextures] = getNGTDMtextures(NGTDM,countValid);
    cellperson = struct2cell(ngtdmTextures);
    len = length(cellperson);
    fprintf(fid, 'NGTDM\n');
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
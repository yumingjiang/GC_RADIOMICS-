function Untitled( dcm_path,nii_path, csv_path )

addpath('C:\Users\95720\Desktop\NIfTI_20140122')
filename = csv_path;
fid = fopen(filename,'a+');

% ***************��ȡ��ʼ��************************************************
%%% dcm�ļ�
I=dicomread(dcm_path); %��ȡͼ��
% I=dicomread('C:\Users\95720\Desktop\ROI ͼ��ͼ\1\ser004img00050.dcm'); %��ȡͼ��
% imagesc(I);%��ʾͼ��
img =double(I);              %���Ҷȼ�ӳ�䵽0~255
low=min(min(img));
high=max(max(img));
maxgray=high-low;           %���㴰��
rate=256/maxgray;
img=I*rate;
img=img+abs(min(min(img)));    %�Ӵ�
dcm_img=uint8(img);           %ת��Ϊ8λ��λͼ���ݸ�ʽ
% imshow(dcm_img)




%%% nii�ļ�
info  = load_nii (nii_path);
% info  = load_nii ('C:\Users\95720\Desktop\ROI ͼ��ͼ\1\1.nii');
nii_img = info.img;
nii_img = imrotate(nii_img, 270,'nearest');
nii_img(:,:)=flipud(nii_img(:,:));
% figure(3);
% imshow(nii_img,[])
% *************************************************************************



% *************ԭͼ��mask����**********************************************
roi_img_w = dcm_img.*uint8(nii_img);
% figure(4)
% imshow(roi_img_w);
% *************************************************************************

% *****************************ԭͼ����ʾroi*******************************
show_img = dcm_img;
[rows , cols ] = size(show_img);
%%% ***********��С��Ͼ���
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

% **************************������roi************************************
rect=[left_line up_line right_line-left_line down_line-up_line];%����ε���С��Ӿ��� 
roi_true=imcrop(dcm_img,rect);%�ü�ROI��Ӿ��Σ��˴�rect=[xmin ymin width hight] 
% imshow(roi_true);
% *************************************************************************


% ******************�󷽲�*************************************************
var_img = double(roi_true);
val_value=var(var_img(:)); 
% *************************************************************************

% *********************��׼��**********************************************
std_img = double(roi_true);
std_value = std2(std_img);
% *************************************************************************

% **************************�����С����ֵ��λ��****************************
max_min_img = roi_true;
max_value = max(max(max_min_img));
[max_row max_col] = find(max_value == max_min_img);
min_value = min(min(max_min_img));
[min_row min_col] = find(min_value == max_min_img);
% *************************************************************************





% *******************************GLCM�������������ء����Ծء����***********
% GLCM = graycomatrix(roi_true,'GrayLimits',[]);
% stats = graycoprops(GLCM,{'energy','contrast','correlation','homogeneity'});

% gray = roi_true
% offsets = [0 1;-1 1;-1 0;-1 -1];
% m = 3; % 3�׻Ҷȼ�
% [GLCMS,SI] = graycomatrix(gray,'GrayLimits',[],'Of',offsets,'NumLevels',m);
% P = GLCMS;
% [kk,ll,mm] = size(P);
% % �Թ��������һ��
% %---------------------------------------------------------
% for n = 1:mm
%    P(:,:,n) = P(:,:,n)/sum(sum(P(:,:,n)));
% end
%  
% %-----------------------------------------------------------
% %�Թ�����������������ء����Ծء����4���������
% %-----------------------------------------------------------
% H = zeros(1,mm);
% I = H;
% Ux = H;      Uy = H;
% deltaX= H;  deltaY = H;
% C =H;
% for n = 1:mm
%    E(n) = sum(sum(P(:,:,n).^2)); %����
%    for i = 1:kk
%       for j = 1:ll
%           if P(i,j,n)~=0
%              H(n) = -P(i,j,n)*log(P(i,j,n))+H(n); %��
%           end
%           I(n) = (i-j)^2*P(i,j,n)+I(n);  %���Ծ�
% 
%           Ux(n) = i*P(i,j,n)+Ux(n); %������Ц�x
%           Uy(n) = j*P(i,j,n)+Uy(n); %������Ц�y
%       end
%    end
% end
% for n = 1:mm
%    for i = 1:kk
%        for j = 1:ll
%             deltaX(n) = (i-Ux(n))^2*P(i,j,n)+deltaX(n); %������Ц�x
%             deltaY(n) = (j-Uy(n))^2*P(i,j,n)+deltaY(n); %������Ц�y
%             C(n) = i*j*P(i,j,n)+C(n);
%        end
%    end
%    C(n) = (C(n)-Ux(n)*Uy(n))/deltaX(n)/deltaY(n); %�����
% end
%  
% %--------------------------------------------------------------------------
% %���������ء����Ծء���صľ�ֵ�ͱ�׼����Ϊ����8ά��������
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
% sprintf('0,45,90,135�����ϵ���������Ϊ�� %f, %f, %f, %f',E(1),E(2),E(3),E(4))  % �������;
% sprintf('0,45,90,135�����ϵ�������Ϊ�� %f, %f, %f, %f',H(1),H(2),H(3),H(4))  % �������;
% sprintf('0,45,90,135�����ϵĹ��Ծ�����Ϊ�� %f, %f, %f, %f',I(1),I(2),I(3),I(4))  % �������;
% sprintf('0,45,90,135�����ϵ����������Ϊ�� %f, %f, %f, %f',C(1),C(2),C(3),C(4))  % �������;
% *************************************************************************
sigma_m = [1.0, 1.5,2.0, 2.5];
fprintf(fid, 'ORI\n');
for count = 1:5
    dcm_img_do = dcm_img;
    if count ~= 1
        %%%% �˲�
        fprintf(fid, '%1.1f�˲�\n',sigma_m(count-1));
        sigma = sigma_m(count-1);
        gausFilter = fspecial('gaussian', [5,5], sigma);
        dcm_img_do = imfilter(dcm_img, gausFilter, 'replicate');
        % figure(2)
        % imshow(dcm_img)
        %%%% 
    end

    [ROIonly,levels] = prepareVolume(dcm_img_do,nii_img,'PETscan',4,3.27,1,5,'Matrix','Lloyd',32);
    
    [GLCM] = getGLCM(ROIonly,levels); %����getGLCM���GCLM����
    [glcmTextures] = getGLCMtextures(GLCM);%����getGCLMtextures�������GCLM����
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
% % ����С���ֽ�
% % % =================================================
% % 
% % ��ͼ�����2���С���ֽ�
% % N = 2 ;
% % [c,s] = wavedec2(ROIonly,N,'db1');
% % 
% % % =================================================
% % ��ȡ�����Ӧ�ĵ�Ƶ�͸�Ƶϵ��
% % % =================================================
% % 
% % ����ԭʼͼ��512��512 ����С���ֽ��1��ά��Ϊ256��256����2��ά��Ϊ128��128
% % ��ȡС���ֽ��е�1���Ƶϵ���͸�Ƶϵ��
% % a_ca1 = appcoef2(c,s,'db1',1);
% % a_ch1 = detcoef2('h',c,s,1);
% % a_cv1 = detcoef2('v',c,s,1);
% % a_cd1 = detcoef2('d',c,s,1);
% % ��ʾ��1��ĸ�����
% % figure(6);
% % subplot(4,4,[3,4,7,8]);imshow(a_ch1,[]);
% % subplot(4,4,[9,10,13,14]);imshow(a_cv1,[]);
% % subplot(4,4,[11,12,15,16]);imshow(a_cd1,[]);
% % 
% % ��ȡ��2��ĵ�Ƶϵ���͸�Ƶϵ��
% % ca2 = appcoef2(c,s,'db1',2);
% % ch2 = detcoef2('h',c,s,2);
% % cv2 = detcoef2('v',c,s,2);
% % cd2 = detcoef2('d',c,s,2);
% % % ��ʾ��2��ĸ�����
% % subplot(4,4,1);imshow(ca2,[]);
% % subplot(4,4,2);imshow(ch2,[]);
% % subplot(4,4,5);imshow(cv2,[]);
% % subplot(4,4,6);imshow(cd2,[]);
% % 
% % 
% % % =================================================
% % �ֱ�Ը�Ƶ�ʳɷֽ����ع�
% % % =================================================
% % 
% % ʹ�õ�2������ع�
% % recon_a1 = wrcoef2('a',c,s,'db1',2);
% % recon_h1 = wrcoef2('h',c,s,'db1',2);
% % recon_v1 = wrcoef2('v',c,s,'db1',2);
% % recon_d1 = wrcoef2('d',c,s,'db1',2);
% % ��ʾ���ع��ĳɷ� ά�ȶ���512��512
% % recon_set = [recon_a1,recon_h1;recon_v1,recon_d1];  
% % figure(7);imshow(recon_set,[]);title('��2��С��ϵ�����ع�');
% % 
% % % =================================================
% % �ع���ԭʼͼ��
% % % =================================================
% % 
% % recon_img = recon_a1+recon_h1+recon_v1+recon_d1;
% % recon_img = mat2gray(recon_img );
% % figure(8);imshow(recon_img );title('�ع�����ԭʼͼ��');
% % WT = wavedec2(ROIonly, 2, 'db4');
% % LL = WT.dec;


fclose(fid);
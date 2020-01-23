clc;
clear;
maindir = 'C:\Users\Yi\Desktop\95720\ROI2';
subdir =  dir(maindir);   % ��ȷ�����ļ���
print(subdir) 
for i = 1 : length( subdir )
    if( isequal( subdir( i ).name, '.' ) || ...
        isequal( subdir( i ).name, '..' ) || ...
        ~subdir( i ).isdir )   % �������Ŀ¼����
        continue;
    end
     
    subdirpath = fullfile( maindir, subdir(i).name, '*.dcm' );
    images = dir( subdirpath );   % ��������ļ������Һ�׺Ϊdcm���ļ�
    dcm_path = fullfile( maindir, subdir(i).name, images(1).name  );
    
    subdirpath = fullfile( maindir, subdir(i).name, '*.nii' );
    images = dir( subdirpath );   % ��������ļ������Һ�׺Ϊnii���ļ�
    nii_path = fullfile( maindir, subdir(i).name, images(1).name  );
    
    fund_path = fullfile( maindir, subdir(i).name);
    feature( dcm_path,nii_path, fund_path );
end
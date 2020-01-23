clc;
clear;
maindir = 'C:\Users\Yi\Desktop\95720\ROI2';
subdir =  dir(maindir);   % 先确定子文件夹
print(subdir) 
for i = 1 : length( subdir )
    if( isequal( subdir( i ).name, '.' ) || ...
        isequal( subdir( i ).name, '..' ) || ...
        ~subdir( i ).isdir )   % 如果不是目录跳过
        continue;
    end
     
    subdirpath = fullfile( maindir, subdir(i).name, '*.dcm' );
    images = dir( subdirpath );   % 在这个子文件夹下找后缀为dcm的文件
    dcm_path = fullfile( maindir, subdir(i).name, images(1).name  );
    
    subdirpath = fullfile( maindir, subdir(i).name, '*.nii' );
    images = dir( subdirpath );   % 在这个子文件夹下找后缀为nii的文件
    nii_path = fullfile( maindir, subdir(i).name, images(1).name  );
    
    fund_path = fullfile( maindir, subdir(i).name);
    feature( dcm_path,nii_path, fund_path );
end
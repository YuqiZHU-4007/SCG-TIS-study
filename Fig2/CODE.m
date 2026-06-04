%% Preprocess resting state data
% Read me
% Edited by Lei Zhao, validated by Jie Chen, supervised by Yiheng Tu
% Add Spm and Conn Toolboxes to path before Running the code
% Steps: 1, Remove first N time points
%        2, Slice timing 
%        3, Realign
%        4, Outliers detection
%        5, Coregister&Segment 
%        6, Normalization 
%        7, Smooth 
%        8, Denoising (detred, regress{'mean signal of White Matter&CSF','friston-24','outliers','Effect of rest'})
%        9, Filter 
%        10, Head motion estimation (Max value of six directions, Power FD at each time points and mean Power FD)
% Inderect normalization is employed in the code (FunImg was coregisted to T1Img then Normalized to MNI space)
%% custom parameters

% Specify your paths and parameters 
%1--- rootpath 
RootPath = 'D:\stim_2\DATA'; % FunImg (named as session 1, session 2 and ......) and T1Img DIR need to be included in the rootpath. put .nii files into img DIRs, 
FunfileName = '4D.nii'; % file name of original funtional.nii file; To avoid multiple .nii files were considered (eg., The residual files of last preprocessing). If you ensure only one .nii file in the dir,you can use *.nii
T1Name = 'T1.nii'; % file name of original T1 .nii file.  To avoid multiple .nii files were considered (eg., The residual files of last preprocessing). If you ensure only one .nii file in the dir,you can use *.nii
disnumber = 10; % discard first N scans before preprocessing

%2--- tamplate and savepath
load('D:\stim_2\LeiIndirect2.mat') % Preprocessed template
SetMat = 'D:\stim_2\DATA\conn_HCP.mat'; % Save your setup in this mat file
RUNPARALLEL = false; % if ture process using paralelle, Default:false
NJOBS = 8; % how many paralells

%3--- The parameters used for fmri preprocessing 
TR = 0.75; % TR
SmoothGaussian = 8; % % smoothing fwhm (mm)
Band = [0.01, 0.08]; % filter Hz

%% discard scans
%Sess001中，会丢弃前N个扫描
cd(RootPath);
[~,SessPaths] = GetDirName(fullfile(RootPath),'Sess*'); % find the number of Session


for i = 1:length(SessPaths)
    
[subnameFun(:,i),subnameFuns(:,i)] = GetDirName(fullfile(RootPath,SessPaths{i})); % get path, each row is a subject, column is a Session
[subnameT1,subnameT1s] = GetDirName(fullfile(RootPath,'T1Img'));  % get path, each row is a subject %其中每一行表示Sub00*的路径

end

NSUBJECTS = size(subnameFun,1); % NSUBJECTS被赋值为subnameFun矩阵的行数，表示Sub的数量。
Nsessions = size(subnameFun,2); % Nsessions被赋值为subnameFun矩阵的列数，表示Sess的数量。

for nsub = 1:NSUBJECTS
    
    for nses = 1:Nsessions
        
        if nses ==1 %检查当前Sess是否为第一个Sess
            [FunctionCeil,myfile] = GetFuncell(subnameFun{nsub,nses}, FunfileName, disnumber);  % only removes scans in the first session
        else
             [FunctionCeil,myfile] = GetFuncell(subnameFun{nsub,nses}, FunfileName, disnumber);% In the remaining sessions, no scans will be removed (disnumber = 0).
        end
        
        matlabbatch{1}.spm.util.cat.vols = FunctionCeil;
        matlabbatch{1}.spm.util.cat.name = ['n',myfile];  % n4D is a new file removing N scans
        matlabbatch{1}.spm.util.cat.dtype = 4;
        matlabbatch{1}.spm.util.cat.RT = NaN;
        spm_jobman('run', matlabbatch); %使用 spm_jobman 函数运行指定的job，其中 matlabbatch 是包含信息
        clear  matlabbatch
    end
end

%% check name homogeneity of T1Img and FunImg
%%通过这段代码，可以确保在每个sess中，每个Sub的功能像和T1匹配（Subname是一致的）。如果不一致，代码将停止执行并显示相应的错误消息。
for i = 1:Nsessions
    
    if ~isequal(subnameFuns(:,i),subnameT1s)

        error(['names of Subject in sess',num2str(i),' and T1Img Must Be Equal']);

    end

end

%% CONN Setup   
%% set preprocessing
clear batch;
batch.filename= fullfile(SetMat);
batch.Setup.isnew= 1;
batch.Setup.nsubjects= NSUBJECTS;
batch.Setup.outputfiles= [0,1,0,0,0,0]; 

% setup condition
batch.Setup.conditions.names={'rest'};                
    for ncond=1
        for nsub=1:NSUBJECTS
            for nses=1:Nsessions             
                batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; 
                batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;
            end
        end
    end     % rest condition (all sessions)

% setup fMRI file (% Point to functional volumes for each subject/session)
% 通过这段代码，为每个sub和sess设置了fMRI的完整路径。这将在后续的预处理步骤中使用。整个过程通过遍历sub和sess，从每个路径中获取以'n'开头的fMRI文件名，并将完整路径存储到batch.Setup.functionals。
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);     

for nsub=1:NSUBJECTS
    
    for nses=1:Nsessions      
            
        thisfile = ls(fullfile(subnameFun{nsub,nses},['n',FunfileName])); % 获取指定路径下以 'n' 开头的功能性 MRI 文件名
        batch.Setup.functionals{nsub}{nses}= fullfile(subnameFun{nsub,nses},thisfile);

    end
    
end 


% setup T1 MRI file
%每个Sub T1完整路径。遍历每个被试，从应路径中获取T1文件名，并将完整路径存储在 batch.Setup.structurals
for nsub = 1:NSUBJECTS
     
    thisfile = ls(fullfile(subnameT1{nsub,1},T1Name));
    STRUCTURAL_FILE{nsub,1} = fullfile(subnameT1{nsub,1},thisfile);
    
end

batch.Setup.structurals= STRUCTURAL_FILE;   

% setup 
batch.Setup.preprocessing.steps= STEPS; % load('F:\ConnBatch\LeiIndirect2.mat')中定义了一系列步骤被存储在变量STEPS中。然后，将STEPS赋值给batch.Setup.preprocessing.steps，以便在CONN中使用这些步骤进行预处理
batch.Setup.RT= TR ;      
batch.Setup.voxelresolution= 1;                          % Volume-based template
batch.Setup.preprocessing.voxelsize_func= 2;          % default 2mm isotropic voxels analysis space if not specify
batch.Setup.preprocessing.fwhm= SmoothGaussian;          % smoothing fwhm (mm)
batch.Setup.preprocessing.art_thresholds= [3, .5]; % position [1]: SD标准差 of Global signal; [2] FD帧间位移

%% CONN Denoising   
%去噪
batch.Denoising.filter= Band;                    % frequency filter (band-pass values, in Hz)
batch.Denoising.detrending= 1;
batch.Denoising.done= 1;                                 % use default denoising step (CompCor, motion regression, scrubbing, detrending)
batch.Denoising.overwrite= 'Yes';
batch.Denoising.confounds.names= {'White Matter','CSF','realignment','scrubbing'}; 
batch.Denoising.confounds.power= {1,1,2,1}; % A parameter '2' adds the quadratic component, because we wanna regress 24-head motion
batch.Denoising.confounds.dimensions= {1,1,inf,inf}; % 5 columns of[1]and[2]include mean white/csf signals and first four components

%% Other 
if RUNPARALLEL
    batch.parallel.N= NJOBS;                             % number of parallel processing batch jobs
end     

batch.Setup.overwrite='Yes';                           
batch.Setup.done=1;

%% measure analysis
% vvAnalysis.done= 1;
% vvAnalysis.overwrite= 1;
% vvAnalysis.measures= {'ALFF';'fALFF'};

%% Run batch
conn_batch(batch);
clear batch;

%% Get motion parameters
mkdir(fullfile(RootPath,'RealignParameter'));%创建一个名为"RealignParameter"的文件夹，用于存储motion parameters文件
for nses = 1:Nsessions %这是一个循环，用于遍历每个Sess
    
    mkdir(fullfile(RootPath,'RealignParameter',SessPaths{nses}));%在"RealignParameter"文件夹下创建名为当前Sess的文件夹，用于存储Sess下的motion parameters文件
    
    for nsub = 1:NSUBJECTS %这是一个嵌套循环，用于遍历每个Sub
        
        outpath = fullfile(RootPath,'RealignParameter',SessPaths{nses},subnameFuns{nsub,nses}); %生成Sub路径
        mkdir(outpath); %创建名为当前Sub的文件夹
        rp = dir(fullfile(subnameFun{nsub,nses},'rp_*')); rp = fullfile(rp.folder,rp.name); copyfile(rp,outpath); %查找名为"rp_*"的motion parameters文件，将其复制到Sub文件夹中。
        rp_value = load(rp); %将motion parameters加载到变量rp_value中，该变量包含了每个时间点的数值
        
        Max_six{nses}(nsub,:) = max(abs(rp_value)); Max_six{nses}(nsub,4:6) = Max_six{nses}(nsub,4:6)*180/pi; %计算每个Sub在当前Session下的motion parameters的最大值，并将角度值转换为度。结果存储在Max_six变量中

        FD_Power{nses,nsub} = GetPower(rp_value);%计算motion power，并将结果存储在FD_Power变量中

        FD_Power_Mean(nses,nsub) = mean(FD_Power{nses,nsub});% 计算motion power的平均值，并将结果存储在FD_Power_Mean变量中。

    end

end
save(fullfile(RootPath,'FD.mat'),'Max_six','FD_Power','FD_Power_Mean'); % 将参数保存到文件中

%% Function
% find subdirs having specific prex
function [files,filename] = GetDirName(path,prex,other)

cd (path);
files = ls ;

if nargin > 1 
    
    files = ls(prex); 
%         files = files(strmatch(prex,files),:)
else 
    
   files(1:2,:)=[];
    
end

files = cellstr(files);
filename = files;
files = cellfun( @(x) fullfile(path, x), files, 'un', 0);
    
if nargin > 2
    
    files = cellfun( @(x) [x, other], files, 'un', 0);
    
end


end

% gather time points from FunImg data
function [FunctionCeil,myfile] = GetFuncell(subnameFun, ImageName, disnumber)
% disnumber: first x time points will be removed
FunctionCeil = cell(1); cd (subnameFun); myfile =  ls(ImageName); 
slicenum = niftiinfo(fullfile(subnameFun,myfile)); slicenum = niftiread(slicenum); TimePoint = size(slicenum); TimePoint = TimePoint(4);

    
    for ii= disnumber+1:TimePoint %time point; 
        frame_order=strcat(',',num2str(ii,'%01d'));
        FunctionCeil(ii-disnumber,:)= GetDirName(subnameFun,ImageName,frame_order);


    end

end

% calculate FD (Power)
function Power = GetPower(rp_value)

RPDiff = diff(rp_value);
RPDiffSphere=[zeros(1,6);RPDiff]; 
RPDiffSphere(:,4:6)=RPDiffSphere(:,4:6)*50; 
Power=sum(abs(RPDiffSphere),2);            

end
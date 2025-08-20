%% PREPROCESSING MOVIE-HINTS

function[] = O1_preprocessing_Movie_HINTS()
spm_dir = 'M:\spm12\'; % complete spm12 folder with backslash at the end
%addpath(genpath(spm_dir));
spm('defaults', 'FMRI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dicom_dir = '...\DICOM\';

study_dir = '...\NIFTI\';
topup_sep_dir = '...\NIFTI';

%SUBJECTS
subjects = ['BHM_437'];

logfile_dir = '...\Logfiles\';
subject_nums = [1];


start_i = 1;
number_subjects = size(subjects, 1);

TR = 1.5;        % enter Repetition Time
number_slices = 63; % enter number of slices
ref_slice = 1;
number_trials = 648;
images2exclude = 2;
fwhm = 6;
HPfilter = 100; % Thats whats recommended.

%for multiband slice_time and reference will be calculated in the slice
%timing step

runNames = {'Scene_4s', 'Scene_12s', 'Scene_36s',...
    'Shot_4s', 'Shot_12s', 'Shot_36s', 'Rest', 'Inverse'};

logCols = {'ID','block','trial_block','trial_total','image','hierarchy','duration',...
    'onset_ms','number_scan','time_scan','img_duration','place','time','location',...
    'arousal_mean','number_people','duration_ratio'};

NiftiToolsDir = '...\NIfTI_20140122\'; %https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
GM_masking_final_images = true; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preprocessing Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% % PREPARTION
make_folders()
dicom2nifti_epis()
dicom2nifti_anat()
rename_niftis()
%%%%%% reorientation (not script-based)

% FIELD MAP CORRECTION
ap_pa_4d()
epis_4d()

%---topup--->FSL Script O2_topup_MH.sh
%---applytopup--->FSL   O2_topup_MH.sh
epis_3d()

% PREPROCESSING
slicetiming()
realignment()
coregistration()
segmentation()
normalisation_epis()
normalisation_anat()
smoothing()

% REGRESSION
get_block_info()
extract_WM_CSF_aCompCor() % according to Nastase & Hasson (2024)
model_specification()
model_estimation()

% AFTER PROCESSING
organize_residuals()
make_GM_mask()
condtions_4D() % for comparison, needs exakt same number of trials in all blocks across conditions
z_scale_epochs() % also applies an GM mask if you want (GM_masking_final_images = true)

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create folders

    function [] = make_folders()
        
        logfiles = dir(logfile_dir);
        logfile_names = {logfiles.name};
        
        for i=start_i:number_subjects
            if ~exist(fullfile(study_dir,subjects(i,:)),'dir')
                fprintf('Creating new subject directory %s.\n',subjects(i,:));
                success = mkdir(fullfile(study_dir,subjects(i,:)));
                if ~success, error('Cannot create new subject directory.\n'); end
            end
            if ~exist(fullfile(study_dir,subjects(i,:),'Anatomy'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Anatomy'));
            end
            if ~exist(fullfile(study_dir,subjects(i,:),'Epis'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Epis'));
            end
            for f = 1:numel(runNames)
                if ~exist(fullfile(study_dir,subjects(i,:),'Epis',runNames{f}),'dir')
                    mkdir(fullfile(study_dir,subjects(i,:),'Epis', runNames{f}));
                end
            end
            
            if ~exist(fullfile(study_dir,subjects(i,:),'Design'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Design'));
            end
            if ~exist(fullfile(study_dir,subjects(i,:),'Model'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Model'));
            end
            
            if ~exist(fullfile(study_dir,subjects(i,:),'Model_13'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Model_13'));
            end
            
            if ~exist(fullfile(study_dir,subjects(i,:),'Behaviour'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'Behaviour'));
            end
            % for topup
            if ~exist(fullfile(study_dir,subjects(i,:),'topup'),'dir')
                mkdir(fullfile(study_dir,subjects(i,:),'topup'));
            end
            if ~exist(fullfile(topup_sep_dir,subjects(i,:),'topup'),'dir')
                mkdir(fullfile(topup_sep_dir,subjects(i,:),'topup'));
            end
            
            %MOVE Outputfile to Behaviour
            logfile_name_i = logfile_names(startsWith(logfile_names, sprintf('%02d', subject_nums(i))));
            logfile_path = fullfile(logfile_dir, logfile_name_i);
            copyfile(logfile_path{1},...
                fullfile(study_dir,subjects(i,:),'Behaviour',...
                [subjects(i,:),'_', logfile_name_i{1}]));
            
        end
        fprintf('Creating Directories:\n Done.\n')
    end

%% Dicom2Nifti Epis
    function [] = dicom2nifti_epis()
        
        foldersDicom= {'BLOCK1', 'BLOCK2', 'BLOCK3',...
            'BLOCK4', 'BLOCK5', 'BLOCK6', 'REST', 'INV'};
        foldersNiftiRest = {'Rest', 'Inverse'};
        
        for i = start_i:number_subjects % for all subjects
            
            disp('Epi Dicom Conversion for: ')
            disp(subjects(i,:)) % print current subject
            
            % Load logfile
            behav_file = fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'*.txt']);
            logfile = dir(behav_file);
            logfile_path = fullfile(logfile.folder, logfile.name);
            
            T = readtable(logfile_path); T.Properties.VariableNames = logCols;
            conditions = T.block(startsWith(T.image, 'trial'));
            order = conditions(1:number_trials:end);
            save(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']),'order');
            
            foldersNifti_i = [order', foldersNiftiRest];
            
            % navigate to correct folder
            allRawFiles = dir(dicom_dir);
            subFolders = {allRawFiles.name}; % all elemts in dicom dir
            subRaw=char(subFolders(startsWith(subFolders,subjects(i,:))));
            
            filesSubRaw = dir([dicom_dir, subRaw]);
            filesSubRaw = {filesSubRaw.name}; % All elements in the found 'subRaw' folder
            expFolder=char(filesSubRaw(startsWith(filesSubRaw, 'BIPSY_MOVIE-HINTS')));
            
            sequSubRaw=dir(fullfile(dicom_dir, subRaw, expFolder));
            sequSubRaw={sequSubRaw.name}; % all elemts in expFodler
            for f = 1:numel(foldersDicom)
                
                rawEpiFolder=char(sequSubRaw(contains(sequSubRaw, foldersDicom{f})));
                
                
                raw_dir_epis = fullfile(dicom_dir, subRaw, expFolder, rawEpiFolder);
                
                % ----------------------------------------------------------- %
                
                epi_dir = fullfile(study_dir, subjects(i,:), 'Epis', foldersNifti_i{f}); % Check spelling
                
                % ----------------------------------------------------------- %
                
                [files_for_conv, ~] = spm_select('FPList',raw_dir_epis); % list all files from raw EPIs folder
                
                % Batch-Editor commands
                jobs{1}.spm.util.import.dicom.data = cellstr(files_for_conv); % load file names in correct format for SPM
                jobs{1}.spm.util.import.dicom.root = 'flat'; % no predefiend folder structure for SPM
                jobs{1}.spm.util.import.dicom.outdir = {epi_dir}; % the specifed output directory
                jobs{1}.spm.util.import.dicom.protfilter = '.*';  % any filter to only use some files?
                jobs{1}.spm.util.import.dicom.convopts.format = 'nii'; % NIFTI format
                jobs{1}.spm.util.import.dicom.convopts.meta = 0; % Should SPM save any metadata of the participants?
                jobs{1}.spm.util.import.dicom.convopts.icedims = 0; % any addtional information?
                
                
                spm_jobman('run', jobs);
                clear jobs % clear jobs after this step for next subject/step
            end
        end
        fprintf('Epi Dicom Conversion:\n Done.\n') % Print, if step is completed
    end

%% Dicom2Nifti Anatomy
% same procedure now for anatomy
    function [] = dicom2nifti_anat()
        
        for i = start_i:number_subjects
            
            disp('Anatomy Dicom Conversion for: ')
            disp(subjects(i,:))
            
            % steps to find the raw anatomy file folder within your folders
            allRawFiles = dir(dicom_dir);
            subFolders = {allRawFiles.name}; % all folders in the folder which contains raw data folder
            subRaw=char(subFolders(startsWith(subFolders,subjects(i,:))));
            % find the raw_data folder in these folders
            filesSubRaw=dir([dicom_dir, subRaw]); % all folders in raw data folder
            filesSubRaw={filesSubRaw.name}; % names of all folders in raw data folder
            expFolder=char(filesSubRaw(startsWith(filesSubRaw,'BIPSY_MOVIE-HINTS'))); % enter beginning of the name of raw data subfolder
            sequSubRaw=dir(fullfile(dicom_dir, subRaw, expFolder));
            sequSubRaw={sequSubRaw.name};
            rawAnatFolder=char(sequSubRaw(startsWith(sequSubRaw,'T1_MPRAGE'))); % enter beginning of the name of mprage raw data subfolder
            % now eith T1_MPRAGE to have the anatomy
            raw_dir_anat = fullfile(dicom_dir, subRaw, expFolder, rawAnatFolder);
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy'); % check spelling
            
            [files_for_conv, ~] = spm_select('FPList',raw_dir_anat);
            
            % Batch commands
            jobs{1}.spm.util.import.dicom.data = cellstr(files_for_conv);
            jobs{1}.spm.util.import.dicom.root = 'flat';
            jobs{1}.spm.util.import.dicom.outdir = {highres_dir};
            jobs{1}.spm.util.import.dicom.protfilter = '.*';
            jobs{1}.spm.util.import.dicom.convopts.format = 'nii';
            jobs{1}.spm.util.import.dicom.convopts.icedims = 0;
            
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Anatomy Dicom Conversion:\n Done.\n')
    end

%% Rename Niftis
% not necessary
    function [] = rename_niftis()
        
        for i = start_i:number_subjects
            
            disp('Rename Niftis for: ')
            disp(subjects(i,:))
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runNames)
                epi_dir = fullfile(epi_root, runNames{f});
                boldData = dir(fullfile(epi_dir,'*.nii'));
                fileNamesBold = {boldData.name}; % list all niftis in Epis folder
                for iFile = 1:numel(fileNamesBold) % looping throgh all of them
                    newName = sprintf(['_bold_HINTS_',runNames{f},'_%05d.nii'],iFile);  % %05d means numebr has 5 digits and will be completed by leading 0 if it has no 5 digits on its own                          % enter favoured names of epi files
                    if ~exist(fullfile(epi_dir,newName))
                        movefile(fullfile(epi_dir,fileNamesBold{iFile}),fullfile(epi_dir,newName));
                        % move the old file to the ne named file, thus creating
                        % the new files
                    end
                end
            end
            
            
            % same procedure for anatomy but only one file so no loop
            % through niftis necessary
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            mprageData = dir(fullfile(highres_dir,'*.nii'));
            fileNameMprage = mprageData.name;
            newName = '_anat_HINTS.nii';                                                   % enter favoured name of anatomy file
            if ~exist(fullfile(highres_dir,newName))
                movefile(fullfile(highres_dir,fileNameMprage),fullfile(highres_dir,newName));
            end
        end
        fprintf('Renaming:\n Done.\n')
    end

%% Reorientation

% that does not work script based
% make sure to also reorient inverted niftis!!

%% the last three epis and inverted epis are saved as one 4d nifti file for FSL

    function [] = ap_pa_4d()
        
        for i = start_i:number_subjects
            
            % Load order
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            lastBlock = runs.order{end};
            
            last_epi_dir = fullfile(study_dir, subjects(i,:), 'Epis', lastBlock);
            
            [epi_filenames_all, ~] = spm_select('FPList',last_epi_dir,'nii');
            epis_for_topup = cellstr(epi_filenames_all(end-2:end,:));
            
            inverse_epi_dir = fullfile(study_dir, subjects(i,:), 'Epis', 'Inverse');
            [inv_filenames_all, ~] = spm_select('FPList',inverse_epi_dir,'nii');
            inv_for_topup = cellstr(inv_filenames_all(end-2:end,:));
            
            %topup_dir = fullfile(study_dir, subjects(i,:), 'topup');
            topup_dir = fullfile(topup_sep_dir, subjects(i,:), 'topup');
            files=[epis_for_topup;inv_for_topup];
            
            jobs{1}.spm.util.cat.vols = cellstr(files);
            jobs{1}.spm.util.cat.name = fullfile(topup_dir, '4Dfortopup.nii');
            jobs{1}.spm.util.cat.dtype = 0;
            jobs{1}.spm.util.cat.RT = 1.5;
            
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('4d fuer topup:\n Done.\n')
    end

%% all functional images are saved as one 4d image for FSL

    function [] = epis_4d()
        
        for i = start_i:number_subjects
            disp(['Making Epis 4D for: ', subjects(i,:)]);
            %topup_dir = fullfile(study_dir, subjects(i,:), 'topup');
            topup_dir = fullfile(topup_sep_dir, subjects(i,:), 'topup');
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runNames)-1
                disp(runNames{f});
                epi_dir = fullfile(epi_root, runNames{f});
                [epis_for_4D, ~] = spm_select('FPList',epi_dir,'^_bold.*\.nii');
                
                jobs{1}.spm.util.cat.vols = cellstr(epis_for_4D);%4d image rauslassen!
                jobs{1}.spm.util.cat.name = [topup_dir, filesep,'epis4topup_', runNames{f},'.nii'];
                jobs{1}.spm.util.cat.dtype = 0;
                jobs{1}.spm.util.cat.RT = 1.5;
                
                spm_jobman('run', jobs);
                clear jobs
            end
        end
        fprintf('epi 4d:\n Done.\n')
    end

%% after topup, convert 4d back to 3d epis

    function [] = epis_3d()
        
        for i = start_i:number_subjects
            %topup_dir = fullfile(study_dir, subjects(i,:), 'topup');
            topup_dir = fullfile(topup_sep_dir, subjects(i,:), 'topup');
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runNames)-1
                epi_dir = fullfile(epi_root, runNames{f});
                epi_zipped = [topup_dir, filesep,'u_epis4topup_', runNames{f},'.nii.gz'];
                gunzip(epi_zipped);
                topuped_4d = epi_zipped(1:(end-3));
                
                jobs{1}.spm.util.split.vol = cellstr(topuped_4d);
                jobs{1}.spm.util.split.outdir = {epi_dir};
                
                spm_jobman('run', jobs);
                clear jobs
            end
        end
        fprintf('epis 4d to 3d:\n Done.\n')
    end

%% Slice Timing

    function [] = slicetiming()
        
        for i = start_i:number_subjects
            
            disp('Slice Timing for: ')
            disp(subjects(i,:))
            
            % define slice_times/slice_order and reference slice from dicom
            % header of any file exept the first
            % specification of folder structure to navigate to teh raw data
            allRawFiles = dir(dicom_dir);
            subFolders = {allRawFiles.name}; % all folders in the folder which contains raw data folder
            subRaw=char(subFolders(startsWith(subFolders,subjects(i,:))));
            % find the raw_data folder in these folders
            filesSubRaw=dir(fullfile(dicom_dir, subRaw)); % all folders in raw data folder
            filesSubRaw={filesSubRaw.name}; % names of all folders in raw data folder
            expFolder=char(filesSubRaw(startsWith(filesSubRaw,'BIPSY_MOVIE-HINTS')));         % enter beginning of the name of raw data subfolder
            % find the correct study folder
            sequSubRaw=dir(fullfile(dicom_dir, subRaw, expFolder)); % all folders in experiment folder
            sequSubRaw={sequSubRaw.name};% names of all folders in experiment folder
            rawEpiFolder=char(sequSubRaw(contains(sequSubRaw, 'BLOCK1')));  % enter beginning of the name of bold raw data subfolder
            % find the one with PP_TR for the functional data
            raw_dir_epis = fullfile(dicom_dir, subRaw, expFolder, rawEpiFolder);
            % define this PP_TR folder as the path for functional raw data
            
            % load in all raw data
            [dicoms_epis, ~] = spm_select('List',raw_dir_epis);
            file_for_order = fullfile(raw_dir_epis, dicoms_epis(15,:));
            
            hdr = spm_dicom_headers(file_for_order);
            slice_times = hdr{1}.Private_0019_1029;
            %            [~,slice_order] = sort(slice_times);
            save(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_slice_times.mat']), 'slice_times');
            %            ref_slice = ceil(median(slice_order));
            ref_time = slice_times(ref_slice);
            %ref_time = slice_times(ref_slice);
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            % reading in all Epis
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            files_for_slicetime = cell(1,numel(runNames)-2);
            for f = 1:numel(runs.order)
                epi_dir = fullfile(epi_root, runs.order{f});
                [filenames, ~] = spm_select('FPList',epi_dir,'^u.*\.nii');
                files_for_slicetime{1,f} = cellstr(filenames((images2exclude+1):end,:));
            end
            
            % Batch commands
            jobs{1}.spm.temporal.st.scans = files_for_slicetime;
            jobs{1}.spm.temporal.st.nslices = number_slices;
            jobs{1}.spm.temporal.st.tr = TR;
            jobs{1}.spm.temporal.st.ta = TR - (TR/number_slices); % enter TA (TR-(TR/slice_number)), sometimes problematic
            % TA is not relevant in slice timing for multiband, since it is
            % already implicitly included in the slice times and thus does
            % not need to be specified
            jobs{1}.spm.temporal.st.so = slice_times;
            jobs{1}.spm.temporal.st.refslice = ref_time;
            
            jobs{1}.spm.temporal.st.prefix = 'a';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Slice Timing:\n Done.\n')
    end

%% Realignment

    function [] = realignment()
        
        for i = start_i:number_subjects
            
            disp('Realignment for: ')
            disp(subjects(i,:))
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            % reading in all Epis
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            files_for_realign = cell(1,numel(runs.order));
            for f = 1:numel(runs.order)
                epi_dir = fullfile(epi_root, runs.order{f});
                [filenames, ~] = spm_select('FPList',epi_dir,'^a.*\.nii');
                files_for_realign{1,f} = cellstr(filenames);
            end
            
            
            % Batch commands
            jobs{1}.spm.spatial.realign.estwrite.data = files_for_realign;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;      % 0.9 is default
            jobs{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % register to first or mean: 1 means mean, 0 means first
            %initially SPM realigns each session to each other by aligning the first volume
            %from each session to the first volume of the first session and subsequently
            %all volumes within each session are aligned to the first volume of that session
            %after that the volumes from the first realignment step are used to create
            %a mean image and than all volumes are aligned to that mean image
            jobs{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            jobs{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % mean image only
            jobs{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            jobs{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            jobs{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            
        end
        fprintf('Realignment:\n Done.\n')
    end


%% Coregistration

    function [] = coregistration()
        
        for i = start_i:number_subjects
            
            disp('Coregistration for: ')
            disp(subjects(i,:))
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            % read in Mean image from Realignment and Anatomy
            % mean image is only in the 1st run Folder ??
            meanepi_dir_1 = fullfile(study_dir, subjects(i,:), 'Epis', runs.order{1});
            find_meanepi = dir(fullfile(meanepi_dir_1,'mean*.nii'));
            meanepi_file = find_meanepi.name;
            
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            find_highres = dir(fullfile(highres_dir, '_anat*.nii'));
            highres_file = find_highres.name;
            
            % define Dar as in Batch
            jobs{1}.spm.spatial.coreg.estwrite.ref = {fullfile(meanepi_dir_1,meanepi_file)}; % reference image: mean epi (where to coreg to)
            jobs{1}.spm.spatial.coreg.estwrite.source = {fullfile(highres_dir,highres_file)}; % source image: anatomy (what to coreg)
            jobs{1}.spm.spatial.coreg.estwrite.other = {''};
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            jobs{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            jobs{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            jobs{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            
        end
        fprintf('Coregistration:\n Done.\n')
    end

%% Segmentation

    function [] = segmentation()
        
        for i = start_i:number_subjects
            
            disp('Segmentation for: ')
            disp(subjects(i,:))
            
            % load anatomy
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            find_highres = dir(fullfile(highres_dir, '_anat*.nii')); % WHY NOT THE COREG?
            highres_file = find_highres.name;
            
            % Batch commands (lots of it is just the parameters for the
            % different materials
            jobs{1}.spm.spatial.preproc.channel.vols = {[highres_dir filesep highres_file]};
            jobs{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            jobs{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            jobs{1}.spm.spatial.preproc.channel.write = [0 1];
            jobs{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm_dir,'tpm','TPM.nii,1')}; % grey matter
            jobs{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(1).native = [1 1]; % FOR DARTEL
            jobs{1}.spm.spatial.preproc.tissue(1).warped = [1 1]; % FOR DARTEL
            jobs{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm_dir,'tpm','TPM.nii,2')}; % white matter
            jobs{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(2).native = [1 1]; % FOR DARTEL
            jobs{1}.spm.spatial.preproc.tissue(2).warped = [1 1]; % FOR DARTEL
            jobs{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm_dir,'tpm','TPM.nii,3')}; % CSF
            jobs{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            jobs{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(3).warped = [1 1]; % save the noramlized, for correction
            jobs{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm_dir,'tpm','TPM.nii,4')}; % bone
            jobs{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            jobs{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm_dir,'tpm','TPM.nii,5')}; % soft tissue
            jobs{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            jobs{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm_dir,'tpm','TPM.nii,6')}; % air/background
            jobs{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            jobs{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            jobs{1}.spm.spatial.preproc.warp.mrf = 1;
            jobs{1}.spm.spatial.preproc.warp.cleanup = 1;
            jobs{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            jobs{1}.spm.spatial.preproc.warp.affreg = 'mni'; % Standard brain
            jobs{1}.spm.spatial.preproc.warp.fwhm = 0;
            jobs{1}.spm.spatial.preproc.warp.samp = 3;
            jobs{1}.spm.spatial.preproc.warp.write = [0 1]; % foward
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Segmentation:\n Done.\n')
    end


%% Normalisation Epis

    function [] = normalisation_epis()
        
        for i = start_i:number_subjects
            
            disp('Writing the normalised Epis for: ')
            disp(subjects(i,:))
            
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            find_deform = dir(fullfile(highres_dir, 'y_*'));
            deform_file = find_deform.name;
            % load Epis with a and deformation field ('y_*')
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runs.order)
                epi_dir = fullfile(epi_root, runs.order{f});
                [files_for_norm , ~] = spm_select('FPList',epi_dir,'^a.*\.nii');
                jobs{1}.spm.spatial.normalise.write.subj.def = {fullfile(highres_dir, deform_file)};
                jobs{1}.spm.spatial.normalise.write.subj.resample = cellstr(files_for_norm);
                jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
                % Voxel size/dimensions
                jobs{1}.spm.spatial.normalise.write.woptions.vox = [2.5 2.5 2.5]; % enter voxel size of normalised epis
                jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
                jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
                % RUN
                spm_jobman('run', jobs);
                clear jobs
                disp(['DONE normalising: ', runs.order{f}]);
            end
        end
        fprintf('Normalisation of Epis:\n Done.\n')
    end

%% Normalisation Anatomy

    function [] = normalisation_anat()
        
        for i = start_i:number_subjects
            
            disp('Writing the normalised Anatomy for: ')
            disp(subjects(i,:))
            
            % Same procedure for Anatomy
            highres_dir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            find_deform = dir(fullfile(highres_dir, 'y_*'));
            deform_file = find_deform.name;
            find_highres = dir(fullfile(highres_dir, 'm_anat*.nii'));
            highres_file = find_highres.name;
            
            jobs{1}.spm.spatial.normalise.write.subj.def = {fullfile(highres_dir, deform_file)};
            jobs{1}.spm.spatial.normalise.write.subj.resample = {fullfile(highres_dir, highres_file)};
            jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
            % now smaller Voxel size
            jobs{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % enter voxel size of normalised anatomy
            jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
            jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Normalisation of Anatomy:\n Done.\n')
    end

%% Smoothing

    function [] = smoothing()
        
        for i = start_i:number_subjects
            
            disp('Smoothing for: ')
            disp(subjects(i,:))
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runs.order)
                epi_dir = fullfile(epi_root, runs.order{f});
                [files_for_smooth, ~] = spm_select('FPList',epi_dir,'^wa.*\.nii');
                
                % Batch commands
                jobs{1}.spm.spatial.smooth.data = cellstr(files_for_smooth);
                jobs{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];  % enter FWHM
                jobs{1}.spm.spatial.smooth.dtype = 0;
                jobs{1}.spm.spatial.smooth.im = 0;
                jobs{1}.spm.spatial.smooth.prefix = 's';
                %RUN
                spm_jobman('run', jobs);
                clear jobs
                disp(['DONE smoothing: ', runs.order{f}]);
            end
        end
        fprintf('Smoothing:\n Done.\n')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract Block Inforamtion from logfiles

    function[] = get_block_info()
        
        for i = start_i:number_subjects
            disp('Get Block Inforamtion for: ')
            disp(subjects(i,:))
            % Load logfile
            behav_file = fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'*.txt']);
            logfile = dir(behav_file);
            logfile_path = fullfile(logfile.folder, logfile.name);
            
            T = readtable(logfile_path);
            T.Properties.VariableNames = logCols;
            
            onsetData = readmatrix(logfile_path, ...
                delimitedTextImportOptions('DataLines',[2,Inf]),'OutputType', 'string'); % read logfile in a form, so onsets are also read
            onsets = onsetData(startsWith(onsetData, 'trigger'));
            trigger_pres = regexp(onsets, '(?<=nr\s)\d+', 'match'); % Numbers after "nr"
            time_pres = regexp(onsets, '(?<=at\s)\d+', 'match'); % Numbers after "at"
            trigger_pres = str2double([trigger_pres{:}]);
            time_pres = str2double([time_pres{:}]);
            onsets = [trigger_pres',time_pres'];
            onsets_blockstart = onsets(1:2:end,:);
            
            subData = T(startsWith(T.image, 'trial'),:);
            order = subData.block(1:number_trials:end);
            
            for b = 1:numel(order)
                subDataBlock = subData(contains(subData.block, order{b}),:);
                subDataBlock.number_scan_corrected = subDataBlock.number_scan - (onsets_blockstart(b,1) -1 + images2exclude);
                disp(subDataBlock.number_scan_corrected(end)-subDataBlock.number_scan_corrected(1)+1);
                subBlocks.(order{b}) = subDataBlock;
            end
            
            save(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_block_info.mat']),'subBlocks');
            
        end
        fprintf('Extraction:\n Done.\n')
    end


%% Extract the top n(= 5)components from the WM and CSF singal
    function [] = extract_WM_CSF_aCompCor()
        % requires fmri denoising by Danilele Mascali
        fmriDenoise_path = 'Q:\fmecklenbrauck\09_WS24_25\Movie-HINTS_Data\fmri_denoising-master';
        addpath(genpath(fmriDenoise_path));
        addpath(NiftiToolsDir);
        name = {'movX','movY','movZ','rotX','rotY','rotZ',...
            'WMcc1','WMcc2','WMcc3','WMcc4','WMcc5',...
            'CSFcc1','CSFcc2','CSFcc3','CSFcc4','CSFcc5'};
        for i = start_i:number_subjects
            disp('Extracting top components from WM and CF signal for: ')
            disp(subjects(i,:))
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            subAnatDir = fullfile(study_dir, subjects(i,:), 'Anatomy');
            
            ROIs = cell(2,1);
            % find ROIs WM & CSF
            for ri = 1:numel(ROIs)
                ROIs{ri,1} = [subAnatDir, filesep, 'wc', num2str(ri+1),'_anat_HINTS.nii'];
            end
            
            % coregister ROIs to EPI space
            clear jobs;
            
            epi_dir = fullfile(epi_root, runNames{1});
            [epis_for_4D, ~] = spm_select('FPList',epi_dir,'^swau.*\.nii');
            
            jobs{1}.spm.spatial.coreg.write.ref = {epis_for_4D(15,:)}; % any swa EPI
            jobs{1}.spm.spatial.coreg.write.source = ROIs;
            jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
            jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
            jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            spm_jobman('run', jobs);
            
            % binarize again in new reference space
            clear jobs
            for ri = 1:numel(ROIs) % all Rois
                jobs{1}.spm.util.imcalc.input = {[subAnatDir, filesep, 'rwc', num2str(ri+1),'_anat_HINTS.nii']};
                jobs{1}.spm.util.imcalc.output = ['brwc', num2str(ri+1),'_anat_HINTS.nii'] ;
                jobs{1}.spm.util.imcalc.outdir = {subAnatDir};
                jobs{1}.spm.util.imcalc.expression = 'i1 > 0.2';
                jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                jobs{1}.spm.util.imcalc.options.dmtx = 0;
                jobs{1}.spm.util.imcalc.options.mask = 0;
                jobs{1}.spm.util.imcalc.options.interp = 1;
                jobs{1}.spm.util.imcalc.options.dtype = 4;
                
                spm_jobman('run', jobs);
                clear jobs;
            end
            
            % get resized, binarized ROIs
            ROIs = cell(1,2); %correct format
            % find ROIs
            for ri = 1:numel(ROIs)
                ROIs{1,ri} = [subAnatDir, filesep, 'brwc', num2str(ri+1),'_anat_HINTS.nii'];
            end
            
            % For all Runs seperately
            for f = 1:numel(runNames)-2
                disp(runNames{f});
                % make files 4D
                epi_dir = fullfile(epi_root, runNames{f});
                [epis_for_4D, ~] = spm_select('FPList',epi_dir,'^swau.*\.nii');
                
                jobs{1}.spm.util.cat.vols = cellstr(epis_for_4D);
                jobs{1}.spm.util.cat.name = [epi_dir, filesep,'swa4D_', runNames{f},'.nii'];
                jobs{1}.spm.util.cat.dtype = 64; % has to be double float (=64)
                jobs{1}.spm.util.cat.RT = 1.5;
                
                spm_jobman('run', jobs);
                clear jobs
                
                % load 4D data
                swa_nii = load_untouch_nii([epi_dir, filesep,'swa4D_', runNames{f},'.nii']);
                
                % extract top 5 aCompCor compontes of both ROIs
                X = fmri_compcor(swa_nii.img, ROIs, [5,5]');
                %Image needs to be 4D, ROIs same 3D dimension,
                % give number components for each ROI
                
                %Load also Movement
                movement_path = dir(fullfile(epi_dir, 'rp*'));
                movement_file = readmatrix(fullfile(epi_dir, movement_path.name));
                
                % combine both into one mat file for easier usage in
                % reegression
                R = [movement_file, X];
                
                save(fullfile(study_dir, subjects(i,:),'Design',...
                    [subjects(i,:),'_', runNames{f}, '_Move_5CompCor_WM_CSF.mat']),...
                    'R');
            end
        end
    end

%% MODEL SPECIFICATION
    function [] = model_specification()
        
        for i = start_i:number_subjects
            tic
            disp('Model specification for:')
            disp(subjects(i,:))
            
            clear jobs;
            %model_dir = fullfile(study_dir, subjects(i,:), 'Model');
            model_dir = fullfile(study_dir, subjects(i,:), 'Model_13');
            
            % define slice_times/slice_order and reference slice from dicom
            sliceTimes = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_slice_times.mat']));
            slice_times = sliceTimes.slice_times;
            [~,slice_order] = sort(slice_times);
            
            
            jobs{1}.spm.stats{1}.fmri_spec.timing.units = 'scans';
            jobs{1}.spm.stats{1}.fmri_spec.timing.RT = TR;
            jobs{1}.spm.stats{1}.fmri_spec.timing.fmri_t = number_slices;
            jobs{1}.spm.stats{1}.fmri_spec.timing.fmri_t0 = find(slice_order == ref_slice);
            jobs{1}.spm.stats{1}.fmri_spec.bases.hrf.derivs = [0 0];            % Canonical HRF, ohne time und dispersion derivative: sonst [1 0] oder [1 1]
            jobs{1}.spm.stats{1}.fmri_spec.cvi = 'AR(1)'; % FAST only for TR < 1.4 better
            jobs{1}.spm.stats{1}.fmri_spec.dir = {model_dir};
            jobs{1}.spm.stats{1}.fmri_spec.global = 'None';
            jobs{1}.spm.stats{1}.fmri_spec.mthresh = 0.8;
            jobs{1}.spm.stats{1}.fmri_spec.mask = {''};
            %jobs{1}.spm.stats{1}.fmri_spec.mask = {%%ADD GM MASK IF NECESSARY};
            jobs{1}.spm.stats{1}.fmri_spec.volt = 1;
            
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            design_dir = fullfile(study_dir, subjects(i,:), 'Design');
            epi_root = fullfile(study_dir, subjects(i,:), 'Epis');
            for f = 1:numel(runs.order)
                epi_dir = fullfile(epi_root,runs.order{f});
                [model_files, ~] = spm_select('FPList', epi_dir, '^swau.*\.nii');  % prefix der vorverarbeiteten EPIS eingeben
                nuisance_reg_file = fullfile(design_dir,  [subjects(i,:),'_', runNames{f}, '_Move_5CompCor_WM_CSF.mat']);
                
                jobs{1}.spm.stats{1}.fmri_spec.sess(f).scans = cellstr(model_files);
                jobs{1}.spm.stats{1}.fmri_spec.sess(f).multi_reg = {nuisance_reg_file};
                jobs{1}.spm.stats{1}.fmri_spec.sess(f).hpf = HPfilter;
            end
            
            spm_jobman('run',jobs)
            
            copyfile(fullfile(model_dir,'SPM.mat'),...
                fullfile(model_dir,'SPM.model'))
            
            clear jobs
        end
        
        disp('Modellspezifikation: DONE ')
        toc
    end


%% Modellschaetzung
%==========================================================================
    function [] = model_estimation()
        for i = start_i:number_subjects
            tic
            disp('Modellschaetzung fuer:')
            disp(subjects(i,:))
            
            model_dir = fullfile(study_dir, subjects(i,:), 'Model_13');
            %model_dir = fullfile(study_dir, subjects(i,:), 'Model');
            jobs{1}.spm.stats{1}.fmri_est.method.Classical = 1;
            jobs{1}.spm.stats{1}.fmri_est.write_residuals = 1;
            jobs{1}.spm.stats{1}.fmri_est.spmmat = {fullfile(model_dir,'SPM.mat')};
            
            spm_jobman('run',jobs)
            copyfile(fullfile(model_dir,'SPM.mat'),...
                fullfile(model_dir,'SPM.estmodel'))
            clear jobs
            
            disp('Modellschaetzung: DONE')
            toc
        end
    end
%% Residuen Umbenennen

    function[] = organize_residuals()
        for i = start_i:number_subjects
            disp('Organization for:')
            disp(subjects(i,:))
            %model_dir = fullfile(study_dir, subjects(i,:),'Model');
            model_dir = fullfile(study_dir, subjects(i,:),'Model_13');
            epis_root = fullfile(study_dir, subjects(i,:),'Epis');
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            
            for f = 1:numel(runs.order)
                if ~exist(fullfile(model_dir,runs.order{f}),'dir')
                    mkdir(fullfile(model_dir,runs.order{f}));
                end
                
                epis_dir = fullfile(epis_root, runs.order{f});
                [files_run_model,~] = spm_select('FPList',epis_dir,'^swau.*\.nii');
                num_runFiles= size(files_run_model,1);
                
                [files_res, ~] = spm_select('List',model_dir,'^Res_.*\.nii');
                
                for iFile = 1:num_runFiles % looping throgh all of them
                    newName = sprintf(['Res_',runs.order{f},'_HINTS_%05d.nii'],(iFile+images2exclude));  % %05d means numebr has 5 digits and will be completed by leading 0 if it has no 5 digits on its own                          % enter favoured names of epi files
                    if ~exist(fullfile(model_dir,runs.order{f},newName),'file')
                        movefile(fullfile(model_dir,files_res(iFile,:)),fullfile(model_dir,runs.order{f},newName));
                        % move the old file to the ne named file, thus creating
                        % the new files
                    end
                end
            end
            
        end
        disp('Organization: DONE')
    end

%% Grey Matter Mask
% ONLY IF ALL SUBJECTS ARE DONE
    function[] = make_GM_mask()
        normGM = cell(number_subjects,1);
        expression_string = '(';
        for i = 1:number_subjects
            normGM{i,1} = fullfile(study_dir, subjects(i,:), 'Anatomy', 'wc1_anat_HINTS.nii');
            expression_string = [expression_string,'i', num2str(i), '+'];
        end
        expression_string = [expression_string(1:(end-1)), ')/',num2str(number_subjects)];
        clear jobs;
        
        jobs{1}.spm.util.imcalc.input = normGM;
        jobs{1}.spm.util.imcalc.output = 'norm_GM_mask.nii';
        jobs{1}.spm.util.imcalc.outdir = {study_dir};
        jobs{1}.spm.util.imcalc.expression = expression_string;
        jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        jobs{1}.spm.util.imcalc.options.dmtx = 0;
        jobs{1}.spm.util.imcalc.options.mask = 0;
        jobs{1}.spm.util.imcalc.options.interp = 1;
        jobs{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % resize to later fit image dimensions
        epi_dir_run1 = fullfile(study_dir, subjects(1,:), 'Epis', runNames{1});
        [func_filenames_wa, ~] = spm_select('FPList',epi_dir_run1,'^wa.*\.nii');
        jobs{1}.spm.spatial.coreg.write.ref = {func_filenames_wa(15,:)}; % 15 ist just any number
        jobs{1}.spm.spatial.coreg.write.source = {fullfile(study_dir,'norm_GM_mask.nii')};
        jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
        jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
        jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'resized_';
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % Smoothing in same dimensions as functional kernel
        jobs{1}.spm.spatial.smooth.data = {fullfile(study_dir, 'resized_norm_GM_mask.nii')};
        jobs{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
        jobs{1}.spm.spatial.smooth.dtype = 0;
        jobs{1}.spm.spatial.smooth.im = 0;
        jobs{1}.spm.spatial.smooth.prefix = ['s',num2str(fwhm),'_'];
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % threshold at 0.25
        jobs{1}.spm.util.imcalc.input = {[study_dir,'/s',num2str(fwhm),'_resized_norm_GM_mask.nii']};
        jobs{1}.spm.util.imcalc.output = 'bin_GM_mask.nii';
        jobs{1}.spm.util.imcalc.outdir = {study_dir};
        jobs{1}.spm.util.imcalc.expression = 'i1 > 0.25'; % based on Nastase et al., 2019
        jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        jobs{1}.spm.util.imcalc.options.dmtx = 0;
        jobs{1}.spm.util.imcalc.options.mask = 0;
        jobs{1}.spm.util.imcalc.options.interp = 1;
        jobs{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run', jobs);
        clear jobs;
        
        fprintf('Grey Matter Mask:\n Done.\n')
    end

%% Concatenate to a single 4D image per condition
    function [] = condtions_4D()
        % thereby leave out 6-10 seconds of data per chunk
        for i = start_i:number_subjects
            disp('4D transformation for:')
            disp(subjects(i,:));
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            % GGF. COPY PART ON LOGFILE FROM extract_WM_CSF_signals
            runsInfo = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_block_Info.mat']));
            
            for f = 1:numel(runs.order)
                disp(runs.order{f});
                curOnset = runsInfo.subBlocks.(runs.order{f}).number_scan_corrected(1);
                curOffset = runsInfo.subBlocks.(runs.order{f}).number_scan_corrected(end);
                
                % leave out first 6-10 seconds
                curOnset_corrected = curOnset + 5;
                % five scans * TR = 7.5 seconds
                
                number_files3D = curOffset - curOnset_corrected +1;
                disp([runs.order{f},': ', num2str(number_files3D)]);
                files3D = cell(number_files3D,1);
                
                %modelRunDir = fullfile(study_dir,subjects(i,:),'Model',runs.order{f});
                modelRunDir = fullfile(study_dir,subjects(i,:),'Model_13',runs.order{f});
                
                m = 1;
                for iFile = curOnset_corrected:curOffset
                    fileName = sprintf(['Res_',runs.order{f},'_HINTS_%05d.nii'],(iFile));
                    pathName = fullfile(modelRunDir, fileName);
                    if exist(pathName, 'file')
                        files3D{m,1} =  pathName;
                    else
                        disp(['Does not exist: ', fileName]);
                    end
                    m = m+1;
                end
                
                clear jobs;
                jobs{1}.spm.util.cat.vols = files3D;
                jobs{1}.spm.util.cat.name = fullfile(modelRunDir, [runs.order{f},'_Res_HINTS.nii']);
                jobs{1}.spm.util.cat.dtype = 64; %FLOAT64 (double precision float value)
                jobs{1}.spm.util.cat.RT = TR;
                
                spm_jobman('RUN', jobs);
                clear jobs;
                disp(['DONE concatenating: ',runs.order{f}]);
            end
        end
        disp('4D transformation: DONE');
    end

%% z_scale_epochs()
% When citing these tools in publications, please reference https://robjellis.net/tools.html​
%z-scored (standardized to zero mean and unit variance) for each voxel and segment

    function [] = z_scale_epochs()
        
        addpath(NiftiToolsDir);
        
        if GM_masking_final_images
            maskPath = fullfile(study_dir, 'bin_GM_mask.nii');
            mask = load_untouch_nii(maskPath);
        end
        
        for i = start_i:number_subjects
            disp('z-scaling Images for:')
            disp(subjects(i,:));
            
            outDir = fullfile(study_dir, subjects(i,:), 'ISC_files_M13');
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            
            runs = load(fullfile(study_dir,subjects(i,:),'Behaviour', [subjects(i,:),'_blockorder.mat']));
            for f = 1:numel(runs.order)
                disp(['RUN: ', runs.order{f}]);
                runDir = fullfile(study_dir, subjects(i,:), 'Model_13', runs.order{f});
                runFile = fullfile(runDir, [runs.order{f},'_Res_HINTS.nii']);%
                
                V = load_untouch_nii(runFile);
                Y = V.img; % Load 4D data (x, y, z, t)
                Vout = V;
                % Get dimensions
                [x_dim, y_dim, z_dim, t_dim] = size(Y);
                
                % Loop through voxels and z-score each voxel’s time series
                b = waitbar(0, 'Starting');
                for x = 1:x_dim
                    for y = 1:y_dim
                        for z = 1:z_dim
                            vox_time_series = squeeze(Y(x, y, z, :));
                            
                            if any(vox_time_series) % Avoid NaNs
                                mVox = mean(vox_time_series);
                                sdVox = std(vox_time_series);
                                if sdVox > 0
                                    Y(x, y, z, :) = (vox_time_series - mVox) / sdVox;
                                else
                                    Y(x, y, z, :) = 0; % If SD is zero, set to zero
                                end
                            end
                        end
                    end
                    waitbar(x/x_dim, b, sprintf('Progress: %d %%', floor(x/x_dim*100)));
                    pause(0.1);
                end
                close(b);
                
                if GM_masking_final_images
                    GM4D = double(repmat(mask.img,[1,1,1, t_dim]));
                    Yfinal = Y.*GM4D;
                else
                    Yfinal = Y;
                end
                
                % Save the z-scored segment
                Vout.img = Yfinal;
                outPath = fullfile(outDir, [runs.order{f},'_z_Res_HINTS.nii']);
                save_untouch_nii(Vout, outPath);
                disp(['DONE: ', runs.order{f}]);
            end
            disp(['DONE subject ', subjects(i,:), ', ', num2str(i),'/',num2str(number_subjects)]);
        end
        disp('z-scaling Images: DONE')
    end

end
%% Extract data from all ISC Maps from all ROIs using MarsBaR

function [] = Extract_ROI_data()

% globals
base_Dir = '...\ROI_analysis'; 
niiTools = '...\NIfTI_20140122'; % https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
addpath(niiTools);
spm_dir = '...\spm12\'; % that includes MarsBar
addpath(genpath(spm_dir));

roiDirs = {'Baldassano_2017_ROIs';'Additional_ROIs'}; 

study_dir = '...\NIFTI';
subjects = ['BHM_437'];

number_subjects = size(subjects,1);

%% SUBFUNCTIONS
tic
make_ROIs_marsy() % the non-spherical ROIs need to be made to MarsBaR roi.mat
extract_data_pairZ_ROIs()
toc
%% ATLAS Creation

    function[] = make_ROIs_marsy()
        
        roiPaths = strcat(base_Dir, filesep, roiDirs);
        
        for p = 1:numel(roiPaths)
            roiDir = roiPaths{p};
            disp(['Writing ROI for: ', roiDirs{p}]);
            roiNiis = dir(fullfile(roiDir, '*.nii'));
            roiNames = {roiNiis.name};
            
            
            for m = 1:numel(roiNames)
                roi = roiNames{m};
                roi = roi(1:(end-4));
                disp(['Making MarsBar RoiMat for: ', roi]);
                imgname = fullfile(roiDir, roiNames{m});
                
                MarsO = maroi_image(struct('vol', spm_vol(imgname), 'binarize',0,'func', 'img'));
                saveroi(MarsO, [roiDir, filesep,roi, '_roi.mat']);
                disp(['DONE ',num2str(m),'/',num2str(numel(roiNames))]);
            end
        end
        
    end
%% MarsBar Extraction

    function[] = extract_data_pairZ_ROIs()
                
        % get all Roi mats
        roiPaths = strcat(base_Dir, filesep, roiDirs);
        allRois = {};
        for p = 1:numel(roiPaths)
            roiCont = dir(fullfile(roiPaths{p}, '*roi.mat')); 
            roiNames = {roiCont.name}';
            roiPath = strcat(roiCont(1).folder, filesep, roiNames); 
            allRois = [allRois; roiPath];
        end
        
        rois = maroi('load_cell', allRois(:));  % make maroi ROI objects
        
        % Step 1: Get filenames with extension (strip folder)
        [~, roiLabel, ~]  = cellfun(@(x) fileparts(x), allRois, 'UniformOutput', false);
        % Step 2: Remove the "_roi" suffix
        roiLabel = erase(roiLabel, '_roi');
        labelTab = cell2table(roiLabel); 
        labelTab.Properties.VariableNames = {'ROI'};
        
        hierarchy = {'Shot','Scene'};
        duration = {'4s','12s','36s'};
        for h = 1:numel(hierarchy)
            hierTab = array2table(repmat(hierarchy(h),numel(rois),1)); 
            hierTab.Properties.VariableNames = {'Hierarchy'};
            for d = 1:numel(duration)
                disp(['Extraction for Condtion: ', hierarchy{h}, duration{d}]);
                durTab = array2table(repmat(duration(d),numel(rois),1)); 
                durTab.Properties.VariableNames = {'Duration'};
                all_subj_nums = 1:number_subjects;
                for i = 1:number_subjects
                    disp(['Exporting for: ', subjects(i,:)]);
                    subjTab = array2table(repmat(string(subjects(i,:)),numel(rois),1)); 
                    subjTab.Properties.VariableNames = {'Subj1'};
                    output = nan(size(rois,1), number_subjects);
                    IscPairDir = fullfile(study_dir,subjects(i,:),'ISC_Pair_maps_M13');
                    imgs =[];
                    for j = 1:number_subjects
                        imgs = [imgs;fullfile(IscPairDir, ...
                            ['ISC_Pair_FishZ_',hierarchy{h},'_',duration{d},'_',subjects(i,:),'_',subjects(j,:),'.nii'])];
                    end
                    non_i = all_subj_nums(all_subj_nums~=i);
                    imgs = imgs(non_i,:);
                    mY = get_marsy(rois{:}, imgs, 'mean');  % extract data into marsy data object
                    y = summary_data(mY); % get summary time course(s)
                    output(:, non_i) = y';
                    save(fullfile(study_dir,subjects(i,:),[subjects(i,:),'_ISCpair_ROIs_',...
                        hierarchy{h},'_',duration{d},'.mat']),'output');
                    outTab = array2table(output);
                    outTab.Properties.VariableNames = cellstr(subjects)';
                    
                    outputTable = [subjTab, hierTab, durTab, labelTab, outTab];

                    writetable(outputTable,...
                        fullfile(study_dir,subjects(i,:),[subjects(i,:),'_ISCpair_ROIs_',...
                        hierarchy{h},'_',duration{d},'.csv']));
  
                    disp(['DONE ',num2str(i),'/',num2str(number_subjects)]);
                end
                disp(['DONE extracting: ', hierarchy{h}, duration{d}]);
            end
        end
    end


end
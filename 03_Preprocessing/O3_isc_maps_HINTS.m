%% ISC map creation
function [] = O3_isc_maps_HINTS()

% This MATLAB script takes all the files in inDir that satisfy a certain
% searchstring then it calculates the average time course and all the
% leave-one-out correlations and saves them in outDir. The script uses the
% NIfTI toolbox, so make sure to install it in your path via:
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
NiftiToolsDir = '...\NIfTI_20140122\';
addpath(NiftiToolsDir);

% Make sure all the 4D files are in the same image space and have the same
% number of time-points or you'll get an error.

% Script written by Christian Keysers
% Adjusted for Movie-HINTS (FM)
% supplemented for pairwise ISC

study_dir = '...\NIFTI\';

subjects = {'BHM_437'};

conditions = {'Scene_4s', 'Scene_12s', 'Scene_36s',...
    'Shot_4s', 'Shot_12s', 'Shot_36s'};% create all paths

number_subjects = numel(subjects);

tic
    isc_pairwise(); %requires fitting brain/GM Mask
toc


%% Pairwise Inter-subject correlation
    function[] = isc_pairwise()
        maskPath = fullfile(study_dir, 'bin_GM_mask.nii');
        mask = load_nii(maskPath);
        mask_voxels = sum(mask.img,'all');
        
        iscPairMapOutDir = fullfile(study_dir, '_ISC_pair_matrices_M13');
        if ~exist(iscPairMapOutDir, 'dir')
            mkdir(iscPairMapOutDir);
        end
        
        % For all Conditions
        for c = 1:numel(conditions)
            % This will setup sum to contain the first image only
            disp(['Calculating ISC maps (Pairwise) for condition: ', conditions{c}]);
            
            iscPairMapOut = NaN(number_subjects,number_subjects,mask_voxels);
            f = waitbar(0,'Please wait...');
            for i=1:number_subjects
                waitbar(i/number_subjects,f,strcat(['starting to process subject ',subjects{i}]));
                fileName = fullfile(study_dir, subjects{i}, 'ISC_files_M13', [conditions{c}, '_z_Res_HINTS.nii']);
                rho_sub=load_nii(fileName);
                [x,y,z,~]=size(rho_sub.img);
                rho = load_nii(fileName,1);
                rho.img=NaN(x,y,z); % maybe need to change to zeros, if there is a mistake
                
                subOut = fullfile(study_dir, subjects{i}, 'ISC_Pair_maps_M13');
                if ~exist(subOut, 'dir')
                    mkdir(subOut);
                end
                
                for j = 1:number_subjects
                    
                    fileName = fullfile(study_dir, subjects{j}, 'ISC_files_M13', [conditions{c}, '_z_Res_HINTS.nii']);
                    tmp=load_nii(fileName);
                    
                    for xi=1:x
                        for yi=1:y
                            for zi=1:z
                                if mask.img(xi,yi,zi)==1 && std(squeeze(rho_sub.img(xi,yi,zi,:)))~=0 && std(squeeze(tmp.img(xi,yi,zi,:)))~=0 % excludes voxels not part of the mask or have no std
                                    rho.img(xi,yi,zi)=corr(squeeze(rho_sub.img(xi,yi,zi,:)),squeeze(tmp.img(xi,yi,zi,:)));
                                end
                            end
                        end
                    end
                    
                    save_nii(rho,strcat([subOut,filesep,'ISC_Pair_corr_',conditions{c},'_',subjects{i},'_',subjects{j},'.nii']));
                    
                    fisherz=rho;
                    fisherz.img=atanh(rho.img);
                    save_nii(fisherz,strcat([subOut,filesep,'ISC_Pair_FishZ_',conditions{c},'_',subjects{i},'_',subjects{j},'.nii']))
                    
                    iscPairMapOut(i,j,:) = fisherz.img(mask.img == 1);
                end
                disp(['DONE: ', subjects{i}]);
            end
            
            iscPairMapOutPath = fullfile(iscPairMapOutDir, ['ISC_Pair_FishZ_',conditions{c},'.mat']);
            save(iscPairMapOutPath,'iscPairMapOut');
            close(f);
        end
    end

end
%% Definition of all 10 bilateral ROIs for Visual Hierachy and supplemental ROI anaylsis
spmDir = '...\spm12';
addpath(genpath(spmDir)); 
addpath('...\NIfTI_20140122') % https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

%% Baldassano et al. (2017) Visual Hierarchy

%% 1.V1 
% Was taken from Anatomy toolbox
% https://github.com/inm7/jubrain-anatomy-toolbox
% Eickhoff S, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K: 
%   A new SPM toolbox for combining probabilistic cytoarchitectonic maps and 
%   functional imaging data. NeuroImage 25(4), 1325-1335, 2005

%% 2. A1: Heschl's Gyrus from Harvard-Oxford
% Dowladed from: https://neurovault.org/collections/262/
%Labels: https://github.com/eglerean/funpsy/blob/master/atlases/HarvardOxford/HarvardOxford-Cortical-Lateralized.xml
cd('...\Heschl_Harvard-Oxford');

gunzip('HarvardOxford-cortl-prob-1mm.nii.gz'); 
label_lat = readstruct('HarvardOxford-Cortical-Lateralized.xml');
label = [label_lat.data.label.Text];
pos = find(contains(label,"Hesch")); 
hemi = {'left','right'};

for i = 1:numel(hemi)
    jobs{1}.spm.util.imcalc.input = {['HarvardOxford-cortl-prob-1mm.nii,',num2str(pos(i))]};
    jobs{1}.spm.util.imcalc.output = ['Heschls_Gyrus_raw_', hemi{i}];
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = 'i1';
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs;
end

hesch_l = niftiread('Heschls_Gyrus_raw_left.nii'); 
hesch_r = niftiread('Heschls_Gyrus_raw_right.nii'); 

max_r = max(max(max(hesch_r)));
max_l = max(max(max(hesch_l)));

maximums = [max_l,max_r];

threshold = 0.25;

for i = 1:numel(hemi)
    jobs{1}.spm.util.imcalc.input = {['HarvardOxford-cortl-prob-1mm.nii,',num2str(pos(i))]};
    jobs{1}.spm.util.imcalc.output = ['Heschls_Gyrus_prob_',hemi{i}];
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = ['(i1./',num2str(maximums(i)), ') > ',num2str(threshold)];
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs;
end


%% 3. hV4: Wang Prob

% dowloaded from: https://napl.scholar.princeton.edu/resources
cd('...\WangProb hV4\ProbAtlas_v4')
gunzip('perc_VTPM_vol_roi7_lh.nii.gz');
gunzip('perc_VTPM_vol_roi7_rh.nii.gz');

hV4_l = niftiread('perc_VTPM_vol_roi7_lh.nii'); 
hV4_r = niftiread('perc_VTPM_vol_roi7_rh.nii'); 

max_l = max(max(max(hV4_l)));
max_r = max(max(max(hV4_r)));

data = {'Left', 'Right';'perc_VTPM_vol_roi7_lh.nii','perc_VTPM_vol_roi7_rh.nii';max_l,max_r};

threshold = 0.25;

for i = 1:size(data,2)
    jobs{1}.spm.util.imcalc.input = data(2,i);
    jobs{1}.spm.util.imcalc.output = ['hV4_025_', data{1,i}];
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = ['(i1./',num2str(data{3,i}), ') > ',num2str(threshold)];
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs;
end

%% 4. Angular gyrus 
% created from Anatomy toolbox: Areas PGa & PGb
% https://github.com/inm7/jubrain-anatomy-toolbox
% Eickhoff S, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K: 
%   A new SPM toolbox for combining probabilistic cytoarchitectonic maps and 
%   functional imaging data. NeuroImage 25(4), 1325-1335, 2005

%% 5. pmDMN
% "POSTERIOR DORSAL DMN" Shirer et al., 2012
% Download from: https://drive.google.com/drive/u/0/folders/1mojjaC3kVbmp7tFU3fWsTQUvuVhMtSXr
cd('...\dorsal_DMN'); 
% posterior medial DMN cluster is cluster number 4
% 1. Cluster is bilateral, we wanted to split to left and rigth hemisphere
% 2. created left and right hemisphere ROI from WFU pick atlas toolbox (Maldjian et al., 2003)
% 3. using SPM's imcalc create binary maps that describe overlap of Shirer
% et al. (2012) pmDMN cluster and left/right hemisphere
pmDMN = fullfile(pwd, '04', '4.nii'); 

%hemispheres to split cluster to left and right are found in the
%_Additional_Files

addFileDir = '...\_AdditionalFiles';
hemi = {'Left.nii';'Right.nii'};


for i = 1:numel(hemi)
    jobs{1}.spm.util.imcalc.input = {pmDMN;fullfile(addFileDir,hemi{i})};
    jobs{1}.spm.util.imcalc.output = ['pmDMN_', hemi{i}];
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = 'i1.*i2';
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs;
end


%% Supplemental ROIs based on areas commonaly realted to the processing of audio-visual narratives aka movies

%%  1. Parahippocampus
% provided from Desikan killiany atlas in "_AddionalFiles"

%% 2. & 3. anterior and posterior hippocampus
cd('.../_AdditionalFiles'); % rpvided von conjunction of HC parcels in Desikan-Killiany atlas

hippo = load_untouch_nii(fullfile(pwd,'Left_Hippocampus_DK.nii')); 
img = hippo.img;% Get voxel indices of ROI (nonzero/binary mask)
voxels = find(img > 0);           % Threshold optional if already binary

% Convert linear indices to subscripts (voxel coordinates)
[x, y, z] = ind2sub(size(img), voxels);

min_y = min(y);
max_y = max(y); 

hippo_lenght = max_y - min_y; 
hipp_middle = floor(min_y + 0.5*hippo_lenght);

img_post = double(img > 0);
img_post(:,hipp_middle:end,:) = 0;
hippo_post = hippo;
hippo_post.img = img_post;

save_untouch_nii(hippo_post, 'Left_post_Hippocampus_DK.nii'); 

img_ant = double(img > 0);
img_ant(:,1:(hipp_middle-1),:) = 0; 
hippo_ant = hippo;
hippo_ant.img = img_ant;
save_untouch_nii(hippo_ant, 'Left_ant_Hippocampus_DK.nii'); 

% Right
hippo = load_untouch_nii(fullfile(pwd,'Right_Hippocampus_DK.nii')); 
img = hippo.img;

% Get voxel indices of ROI (nonzero/binary mask)
voxels = find(img > 0);           % Threshold optional if already binary

% Convert linear indices to subscripts (voxel coordinates)
[x, y, z] = ind2sub(size(img), voxels);

min_y = min(y);
max_y = max(y); 

hippo_lenght = max_y - min_y; 
hipp_middle = floor(min_y + 0.5*hippo_lenght);

img_post = double(img > 0);
img_post(:,hipp_middle:end,:) = 0;
hippo_post = hippo;
hippo_post.img = img_post;

save_untouch_nii(hippo_post, 'Right_post_Hippocampus_DK.nii'); 

img_ant = double(img > 0);
img_ant(:,1:(hipp_middle-1),:) = 0; 
hippo_ant = hippo;
hippo_ant.img = img_ant;
save_untouch_nii(hippo_ant, 'Right_ant_Hippocampus_DK.nii'); 

% Extract transformation components
pixdim = hippo.hdr.dime.pixdim(2:4);     % voxel size
qoffset = [hippo.hdr.hist.qoffset_x, ...
           hippo.hdr.hist.qoffset_y, ...
           hippo.hdr.hist.qoffset_z];    % origin offset

% Build the affine matrix (assuming no rotation/skew — standard RAS orientation)
T = [-pixdim(1), 0,           0,           qoffset(1);  % flip x if needed
     0,         pixdim(2),   0,           qoffset(2);
     0,         0,           pixdim(3),   qoffset(3);
     0,         0,           0,           1];


T_inv = inv(T);
mni_coord = [0; -21; 0; 1];     % MNI (x, y, z, 1)
voxel_coord = T_inv * mni_coord;
voxel_coord = round(voxel_coord);  % round to voxel index

voxel_thresh_y = voxel_coord(2); 

coronal_gap = 2; 

% Left
hippo = load_untouch_nii(fullfile(pwd,'Left_Hippocampus_DK.nii')); 
img = hippo.img;% Get voxel indices of ROI (nonzero/binary mask)

img_post = double(img > 0);
img_post(:,(voxel_thresh_y+coronal_gap):end,:) = 0;
hippo_post = hippo;
hippo_post.img = img_post;

save_untouch_nii(hippo_post, 'Left_post-21gap_Hippocampus_DK.nii'); 

img_ant = double(img > 0);
img_ant(:,1:(voxel_thresh_y-coronal_gap),:) = 0; 
hippo_ant = hippo;
hippo_ant.img = img_ant;
save_untouch_nii(hippo_ant, 'Left_ant-21gap_Hippocampus_DK.nii'); 

% Right
hippo = load_untouch_nii(fullfile(pwd,'Right_Hippocampus_DK.nii')); 
img = hippo.img;

img_post = double(img > 0);
img_post(:,(voxel_thresh_y+coronal_gap):end,:) = 0;
hippo_post = hippo;
hippo_post.img = img_post;

save_untouch_nii(hippo_post, 'Right_post-21gap_Hippocampus_DK.nii'); 

img_ant = double(img > 0);
img_ant(:,1:(voxel_thresh_y-coronal_gap),:) = 0; 
hippo_ant = hippo;
hippo_ant.img = img_ant;
save_untouch_nii(hippo_ant, 'Right_ant-21gap_Hippocampus_DK.nii'); 


%% 4. Temporo-parietal junction
% Connectivity based parcellation atlas (Mars et al., 2012)
%Downloaded from: http://www.rbmars.dds.nl/CBPatlases.htm
cd('...\TPJ Mars_2012\MarsTPJParcellation');

gunzip('TPJ_thr25_1mm.nii.gz');

 jobs{1}.spm.util.imcalc.input = {'TPJ_thr25_1mm.nii,1'};
    jobs{1}.spm.util.imcalc.output = 'TPJ_thr25_Mars_Ant_Right.nii';
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = 'i1 > 0';
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs;
   
 jobs{1}.spm.util.imcalc.input = {'TPJ_thr25_1mm.nii,2'};
    jobs{1}.spm.util.imcalc.output = 'TPJ_thr25_Mars_Post_Right.nii';
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = 'i1 > 0';
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs
    
 jobs{1}.spm.util.imcalc.input = {'TPJ_thr25_1mm.nii,1';'TPJ_thr25_1mm.nii,2'};
    jobs{1}.spm.util.imcalc.output = 'TPJ_thr25_Mars_Right.nii';
    jobs{1}.spm.util.imcalc.outdir = {pwd};
    jobs{1}.spm.util.imcalc.expression = '(i1 + i2) > 0';
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 64;
    
    spm_jobman('run', jobs); 
    clear jobs

% only in right Hemisphere
% mirror
tpj = load_untouch_nii(fullfile(pwd, 'TPJ_thr25_Mars_Right.nii')); 
img = tpj.img;
imagesc(img(:,:,100)); 

img_flipped = flip(img,1); 
imagesc(img_flipped(:,:,100)); 

img = double(img > 0); 
img_flipped= double(img_flipped > 0); 

tpj_right = tpj; tpj_left = tpj;

tpj_right.img = img;
tpj_left.img = img_flipped;

save_untouch_nii(tpj_right, 'Right_TPJ_Mars.nii'); 
save_untouch_nii(tpj_left, 'Left_TPJ_Mars.nii'); 


%% 5. medial preforntal cortex

% cut at x = +-15 (Córcoles-Parada et al., 2017)
cd('.../_Additional_Files');

% Left
sfg = load_untouch_nii(fullfile(pwd,'Left_SFG_lasuanne250D1_3_4_5.nii')); 
img = sfg.img;% Get voxel indices of ROI (nonzero/binary mask)
voxels = find(img > 0);           % Threshold optional if already binary

% Convert linear indices to subscripts (voxel coordinates)
[x, y, z] = ind2sub(size(img), voxels);

min_x = min(x);
max_x = max(x); 

% Extract transformation components
pixdim = sfg.hdr.dime.pixdim(2:4);     % voxel size
qoffset = [sfg.hdr.hist.qoffset_x, ...
           sfg.hdr.hist.qoffset_y, ...
           sfg.hdr.hist.qoffset_z];    % origin offset

% Build the affine matrix (assuming no rotation/skew — standard RAS orientation)
T = [pixdim(1), 0,           0,           qoffset(1);  % flip x if needed
     0,         pixdim(2),   0,           qoffset(2);
     0,         0,           pixdim(3),   qoffset(3);
     0,         0,           0,           1];


T_inv = inv(T);
mni_coord = [-15; 0; 0; 1];     % MNI (x, y, z, 1)
voxel_coord = T_inv * mni_coord;
voxel_coord = round(voxel_coord);  % round to voxel index

voxel_thresh_x = voxel_coord(1); 

img_cut = double(img > 0);
img_cut(1:(voxel_thresh_x-1),:,:) = 0;
sfg_l_cut = sfg;
sfg_l_cut.img = img_cut;

save_untouch_nii(sfg_l_cut, 'Left_mPFC_lausanne250D1_3_4_5.nii');

% Right
sfg = load_untouch_nii(fullfile(pwd,'Right_SFG_lasuanne250D1_2_3_4.nii')); 
img = sfg.img;% Get voxel indices of ROI (nonzero/binary mask)
voxels = find(img > 0);           % Threshold optional if already binary

% Convert linear indices to subscripts (voxel coordinates)
[x, y, z] = ind2sub(size(img), voxels);

min_x = min(x);
max_x = max(x); 

% Extract transformation components
pixdim = sfg.hdr.dime.pixdim(2:4);     % voxel size
qoffset = [sfg.hdr.hist.qoffset_x, ...
           sfg.hdr.hist.qoffset_y, ...
           sfg.hdr.hist.qoffset_z];    % origin offset

% Build the affine matrix (assuming no rotation/skew — standard RAS orientation)
T = [pixdim(1), 0,           0,           qoffset(1);  % flip x if needed
     0,         pixdim(2),   0,           qoffset(2);
     0,         0,           pixdim(3),   qoffset(3);
     0,         0,           0,           1];


T_inv = inv(T);
mni_coord = [15; 0; 0; 1];     % MNI (x, y, z, 1)
voxel_coord = T_inv * mni_coord;
voxel_coord = round(voxel_coord);  % round to voxel index

voxel_thresh_x = voxel_coord(1); 

img_cut = double(img > 0);
img_cut(voxel_thresh_x:end,:,:) = 0;
sfg_r_cut = sfg;
sfg_r_cut.img = img_cut;

save_untouch_nii(sfg_r_cut, 'Right_mPFC_lausanne250D1_2_3_4.nii');




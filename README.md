# Time Points or Plot Points?

This repository contains all the code used in data collection, preparation and analysis in the paper *Time Points or Plot Points - Are movies processed according to their temporal duration or their underlying content structure?* (Mecklenbrauck \& Schubotz, in preparation).

Requirements:

* MATLAB, R(Studio)
* AFNI, FSL, SPM12
* [Tools for NifTI and ANALYZE Image](https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) by Jimmy Shen
* [WFU pickatlas Toolbox](https://www.nitrc.org/projects/wfu_pickatlas/) for SPM12 (if you are interested in a MNI version of the Desikan-Killiany atlas incl. the lausanne-sub-parcellations ([Cammoun et al., 2012](https://doi.org/10.1016/j.jneumeth.2011.09.031))for the WFU toolbox feel free to contact me (f\_meck01@uni-muenster.de)).
* [MarsBaR Toolbox](https://marsbar-toolbox.github.io/download.html) for SPM12
* [Anatomy Toolbox](https://github.com/inm7/jubrain-anatomy-toolbox)

If you want to download a minimally preprocessed example participant (ran O1\_preprocessing\_Movie\_HINTS.m (Part 1)), all pairwise-ISC Maps, or the data for the ROI analyis, click [here](https://uni-muenster.sciebo.de/s/TLFLHTP4ysQgayr).

## Contents

### 01\. Stimuli

The stimuli in this study were created from frames of the movie *Forrest Gump* (Zameckis, 1994). We used the publically available annotation from the [StudyForrest project](https://www.studyforrest.org/about.html) for the selection and assignment of all scenes and shots to the conditions of the experiment. The code describes our selection process from movie to extracted frames.

#### 01\_Scene\_Shot\_Selection.R

Requirements from StudyForrest project: Cuts, Depicted Locations, and Temporal Progression ([Häusler \& Hanke, 2016](https://f1000research.com/articles/4-92/v1)) and Portrayed Emotions ([Labs et al., 2015](https://f1000research.com/articles/4-92/v1)) as well as all files from the AdditionalFiles directory (will be referenced in the Code)

This code loads all necessary data, adjust it so it matches the movie perfectly and then combines scene/shot information with the information of the emotional arousal rating, so each scene and shot has a description of the people depicted, how many people and what there emotional valence is within this movie segment.
Then it selects only shots that that show at least one person, are only a single camera shot or ar part of scenes that were excluded. Scenes, accordingly, contain also at least one person but are also at least made up by at least 2 shots. Additionally, we excluded extremely short or long scenes.
Then the Euclidian distance of emotional arousal, duration and people per scene was used to match similar scenes and shots. Then the scenes and shots were distributed between the 3 duration conditions, minimizing differences in the distributions of arousal, duration and people per scene as well as controlling for similar amounts of Indoor or Outdoor and Day or Night segments. This matching and distribution procedure was repeated 100 times. We then selected the option with the least amount of differences between scenes and shots in perceived tempo when applying the Duration manipulation.

#### 02\_Balacing\_Blocks\_Order.R

The selected scenes and shots were then ordered per condition controlling for average forward and backward jumps through the movie, that two segments were not back-to-back in the original movie and that adjacent segments were from different locations.

#### 03\_automaticCutting.py

Using [ffmpeg](https://ffmpeg.org/) we used the HH:MM:SS:FF formatted time codes of the selected segments to automatically cut clips from the movie.
**Attention**: StudyForrest used a cut of the movie that excluded scenes with discriminating scenes, we had to manually create this cut, however du to copyrights, we are not able to share this version on GitHub.

#### O4\_extract\_isochronic\_frames.m

Loading a clip and based on the duration condition (4s, 12s or 36s) isochonically exporting 12, 36 or 108 frames from each clip.

### 02\. Experiment

This directory contain both the code for the cover task as well as all code to run the experiment in Neurobs Presentation.

#### O1\_Block\_orders.m

Unfortunately, a little manual as it has to be rerun, until enough block orders for 48 participants were created. Creates 6 orders that start with 1:6 and have every number 1:6 once in every position. Does that multiple times and selects 8 orders that are all different and still balance the position of every block in the order. Checks balancing and save the final order

#### O2\_make\_block\_quest\_orders.m

For each subject, this script finds a balanced order of the 6 trials in each block depicts mirrored frames. Thereby, controlling that each block has equally often 2, 3 or 4 trials with mirrored images, each participants sees two locks with 2, 3, or 4 mirrored frames, and each individual frame is equally often the one to be flipped across participants.

#### Presentation experiment

Presentation experiment (kann ich gerade nicht testen haben keine Presentation Lizenz)

### 03\. Preprocessing

#### O1\_preprocessing\_Movie\_HINTS.m (Part 1)

Part 1 describes the preparation of the raw DICOM files in SPM up to the point we then field map corrected the data using topup in FSL.

1. PREPARTION:

   * make\_folders(): creates folder structure
   * dicom2nifti\_epis() and dicom2nifti\_anat(): converts DICOM to NIfTI files
   * rename\_niftis(): relabels nifty files, because SPM filenames are too long, can cause issues
   * Manual reorientation: checking data and rigidly reorienting  anatomical and functional files to MNI template

2. FIELD MAP CORRECTION (also still preparation)

   * ap\_pa\_4d(): creates 4D epi with 6 images (last 3 epis from scan and last 3/5 reverse phase encoding direction images)
   * epis\_4d(): makes all blocks separate 4D epi files to be handled by FSL

#### O2\_topup\_MH.sh

We then switch to FSL for topup field map correction. The code was run on a Linux machine. It estimates the field map correction from AP-PA-phase encoding direction file und then applies correction to all blocks separately.

#### O1\_preprocessing\_Movie\_HINTS.m (Part 2)

Part 2 then continues with the preprocessing and data preparation in SPM.

2. FIELD MAP CORRECTION (also still preparation)

   * epis\_3d(): makes topupped epis to 3D files again to be handled by SPM

3. PREPROCESSING

   * slicetiming(): W used first slice as reference, since due to multislice acquisition no other slice stood out, so we left it at default.
   * realignment(): For multiple blocks/runs SPM first aligns each run’s first volume to the first volume of the first run, then all images are aligned to the first image of their respective block, then a mean image is created and all images are aligned to that mean image. The mean image necessary for coregistration.
   * coregistration(): Aligning anatomy to mean image of functional data.
   * segmentation(): Segmenting the coregistered anatomy into grey matter, white matter, CSF, bone and skin. Grey matter, white matter and CSF are normalized as well to MNI space using estimated deformation field as they will be used for data extaction necessary for then uisance regression.
   * normalisation\_epis(): Normalizing all the function files to MNI space and resampling them to 2.5 x 2.5 x 2.5 mm³ voxels.
   * normalisation\_anat(): Normalizing bias field corrected anatomy (not that necessary)
   * smoothing(): Smoothing with recommended the 6 x 6 x 6 mm³ Gaussian smoothing kernel.

4. REGRESSION

   * get\_block\_info(): Reading log-files from experiment to get trigger information to link experiment events to data. Saving block order for each participant
   * extract\_WM\_CSF\_aCompCor(): Coregistering normalized WM and CSF maps to EPI space (is smaller), binarizing the maps again and making preprocessed functional data of each block/run to one 4D file. fmir\_compcor is then used to extract first 5 components of both the WM and CSF signal across the run time per block. Then this function combines data with movement data of that block and create one nuisance regressors file (“…\_Move\_5CompCor\_WM\_CSF.mat”) per block. This models follow the recommendation of [Nastase \& Hasson (2024)](https://github.com/snastase/isc-confounds).
   * model\_specification(): Using nuisance regressor in noise only model is SPM and applying a 1/100 Hz high pass filter.
   * model\_estimation(): Estiamting the model and importantly saving the residuals.

5. AFTER PROCESSING

   * organize\_residuals(): Sorting the residuals to each block.
   * make\_GM\_mask(): Creating a GM mask from the average normalized grey matter across participants and coregistering it to the size of the residuals.
   * condtions\_4D(): Creating single 4D NIfTI file from the residuals for each block/condtion, while additionally excluding the first 5 images per block. Make sure every subject has the exact same number of time points/volumes.
   * z\_scale\_epochs(): Voxelwise normalization of data. Function also applying gray matter mask to data (if activated). The saves processed block-wise 4D files ready for ISC anaylsis.

#### O3\_isc\_maps\_HINTS.m

1. Code by Christina Keysers adapted from loo-isc to pairwise-isc
2. Requires: Gray Matter mask, [NiftiTools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
3. Loads prepared 4D image of a subject, correlates it voxelwise (only gray matter mask voxel) with every subject
4. Saves correlation file for each pair plus a mat-file for all voxel all correlations (N\_subject x N\_Subject x N\_GMvoxels)

Issues: Code also calculates maps for subject with itself and runs very slowly…

### 04\. Analysis

The Analysis uses the 3dISC framework in AFNI ([Cheng et al., 2017](http://dx.doi.org/10.1016/j.neuroimage.2016.08.029)), so the code is created to run on a Linux machine.

#### 01\_AFNI\_contrasts.sh

Using 3dcalc to separately calculate all contrasts while simultaneously applying Fisher’s z transformation and averaging across conditions if applicable.

#### 02\_MakeDataTables.R

To run 3dISC AFNI require specifically formatted data table files. These files were created for the condition wise averages as well as the contrasts.

#### 03\_\[1-3]\_Runc3dISC…

Separate script for running the 3dISC model estimation based on the data table files using a very simple intercept only multilevel model with crossed random effects for the two subjects involved in each pair:
$z\_{AB} \\sim 1 + (1|\\text{SubjA}) + (1|\\text{SubjB})$

#### 04\_1\_calcCluster\_Conditions.sh

1. Estimating smoothness using 3dFWHMx
2. Simulating cluster size using estimated smoothness (acf) at NN = 3, cluster forming p-threshold = 0.001 and FWE correct cluster size alpha = 0.05 using a bisided test direction  3dClustSim
3. Lastly, using 3dClusterize und threshold the average condition maps with the estimated cluster thresholds

#### 04\_\[2/3]*calcCluster*...

4. Additional steps after the procedure in 04\_1\_calcCluster\_Conditions.sh
5. Using whereami to label all cluster maxima
6. Custom function to determine local maxima per cluster

   * Taking the clusterized map
   * Creating mask from each cluster
   * Using 3dExtrema to look for (5) extreme t-values in each mask that are at least 8 voxels apart (inspired by SPM local maxima definition)

#### O5\_MakeOutputDataTable.m

Based on the text reports for cluster and local maxima
Creating excel output table that sort the local maxima to the cluster maxima and writes out voxel number, location (x, y, z in MNI space) of centre of the cluster, extension of cluster in x, y and z direction, mean value in cluster, SD in cluster, maximal value in cluster as well as local maxima and location of cluster and local maxima (x, y, z in MNI space).

#### ROI analysis

##### O1\_ROI\_defintions.m

Description and code for the creation of all 20 lateralized ROIs across visual hierarchy and supplemental ROIs.
Download links to sources (also provided in Code):

* [A1/Heschl’s Gyrus](https://neurovault.org/collections/262/) plus [Label](https://github.com/eglerean/funpsy/blob/master/atlases/HarvardOxford/HarvardOxford-Cortical-Lateralized.xml)
* [Area hV4](https://napl.scholar.princeton.edu/resources)
* [posterior medial default mode network (pmDMN)](https://drive.google.com/drive/u/0/folders/1mojjaC3kVbmp7tFU3fWsTQUvuVhMtSXr)
* [temporo-parietal junction](http://www.rbmars.dds.nl/CBPatlases.htm)

Rest was either extracted from Desikan-Killiany atlas or the [Anatomy Toolbox](https://github.com/inm7/jubrain-anatomy-toolbox).

##### O2\_Extract\_ROI\_data.m

Using MarsBar to extract the average value from all ROIs from all Fisher’s z-scaled pairwise-ISC maps of each subject pair in each condition.

##### O3\_isc\_Roi\_analysis\_plotting.R

1. Loading all extracted ROI data
2. Modelling in each ROI:
   $\\bar{z}\_{AB,ROI} \\sim \\text{Content Level} \\times \\text{Duration} + (1|\\text{SubjA}) + (1|\\text{SubjB})$

   * Mirroring 3dISC with crossed random effects

3. Visual inspection of Model fit
4. Parametric bootstrapping and saving FDR corrected post-hoc tests
5. Plotting based on bootstrapped values

### 05\. Graphics

This directory specifically refers to the result figures of the whole-brain analysis (Figure 3 – 5). We applied the “Highlight, don’t Hide” principle ([Taylor et al., 2023](https://doi.org/10.1016/j.neuroimage.2023.120138)). Thus, clusters with non-significant p-values are not omitted but displayed at lower opacity, while significant cluster are depicted at full opacity with a black border. While the paper presents the code in volume space we created a version in surface space.

#### HighlightSurfaceMapVisualizer.py

This is the python function that has to be imported.
The code requires:

* the unthresholded results map depicted the raw contrast values
* a map with (t-)statistic values necessary for evaluating opacity by “significance”
* thresholded contrast image, so the program knows which cluster are actual significant

It than determines color and through the raw contrast values, opacity by the absolute statistical values and significant cluster from the thresholded map, map all maps on one surface. Additionally, the program provides the option to select an inflated or non-inflated (inflated = \[True,False]) surface of different resolutions as well as the view and hemisphere you want to have plotted. Default is all views of the selected hemisphere.

Library import requirements:

* numpy
* os
* warnings
* nilearn (surface, datasets)
* pyvista
* matplotlib

#### 01\_runHSMV\_bothHemi.py

Example how the graphics for the paper were created.
For all contrasts of interest select the respective input files (correlation-map, t-statistics-map, significant-cluster-map). The iterate through both hemispheres and the inflation True or False options. We used the surface template ‘fsaverage7’ and the afni\_hotcold colormap that is included as default in HighlightSurfaceMapVisualizer.py. The selection views = None enables that all 6 views (front, back, left, right, top, bottom) are plotted.


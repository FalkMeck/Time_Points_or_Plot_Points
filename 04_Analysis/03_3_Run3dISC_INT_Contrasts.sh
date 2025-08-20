#!/bin/bash

studyDir=/.../NIFTI


contrasts=("INT01_12sSceneShot_04sSceneShot" "INT02_36sSceneShot_04sSceneShot" "INT03_36sSceneShot_12sSceneShot") 


for con in "${contrasts[@]}";do
	
	echo ${con}
	
	cd ${studyDir}/_AFNI_Analysis/${con}

	dataTable=${studyDir}/_AFNI_Analysis/${con}/isc_datatable_${con}.txt
	#echo ${dataTable}

	prefix=Con_${con}

	3dISC -prefix $prefix -jobs 8 \
		 -mask ${studyDir}/bin_GM_mask.nii \
		 -model '1+(1|Subj1)+(1|Subj2)' \
		 -dataTable @${dataTable}
		 
	3dAFNItoNIFTI -prefix ${prefix}_corr.nii ${prefix}+orig'[0]'
	3dAFNItoNIFTI -prefix ${prefix}_tstat.nii ${prefix}+orig'[1]'


done

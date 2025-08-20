#!/bin/bash

studyDir=/.../NIFTI


contrasts=("01_Shot_Scene" "02_04s_12s" "03_04s_36s" "04_12s_36s")


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

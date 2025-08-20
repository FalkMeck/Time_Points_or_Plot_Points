#!/bin/bash
set -e

# Run3dIsc_MH_v02.sh
#
# Created 2025-04-15 by FM.

studyDir=/media/f_meck01/Movie_HINTS/Movie_HINTS_Analysis/NIFTI


outDir=${studyDir}/ISC3D_Con
if [ ! -d ${outDir} ]; then
  mkdir ${outDir}
fi

cd ${outDir}


#echo @{$dataTable}
# Other possible options

hierarchy=("Shot" "Scene")
duration=("4s" "12s" "36s")



for h in "${hierarchy[@]}";do
for d in "${duration[@]}";do


dataTable=${studyDir}/isc_con_table_${h}_${d}_2025_04_22.txt
echo ${dataTable}

prefix=${h}_${d}


3dISC -prefix $prefix -jobs 8 \
	 -mask ${studyDir}/bin_GM_mask.nii \
	 -model '1+(1|Subj1)+(1|Subj2)' \
	 -dataTable @${dataTable}
done
done

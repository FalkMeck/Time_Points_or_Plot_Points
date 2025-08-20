#!/bin/bash
set -e

# Created 2025-04-23 by FM.


# Variables
studyDir=/.../NIFTI

subjects=("BHM_437")

hierarchy=("Shot" "Scene")
duration=("4s" "12s" "36s")

contrasts=("INT01_12sSceneShot_04sSceneShot" "INT02_36sSceneShot_04sSceneShot") #"INT03_36sSceneShot_12sSceneShot")
contrast_indices=("0 1 3 4" "0 2 3 5" "1 2 4 5")

maskFile=${studyDir}/bin_GM_mask.nii

outRoot=${studyDir}/_AFNI_Analysis/ISC3D_Contrasts
mkdir -p ${outRoot}

#Clusterize values
neighbors=3
pThresh=0.001
#alpha=0.05
testDirect="bisided" #ATTENTION: If changed change in clusterize_map aswell

estimate_smoothness() {

for h in "${hierarchy[@]}";do
	for d in "${duration[@]}";do

		outDir=${outRoot}/${h}_${d}_FWHM
		mkdir -p ${outDir}


		cd ${outDir}

		for subject in "${subjects[@]}";do

			in_file=${studyDir}/${subject}/ISC_files/${h}_${d}_z_Res_HINTS.nii

			      #Check input file exists
			      if [ ! -f "${in_file}" ]; then
				echo "** WARNING: Missing file: ${in_file}"
				continue
			      fi

			echo ${in_file}

			3dFWHMx -acf -mask ${maskFile} \
				-input ${in_file} > ${subject}_acf.1D
				#-demed -unif \
				#-out ${subject}_acf.1D \


			sed -n '2p' ${subject}_acf.1D >> group_ACF_ests.1D

		done
	done
done
}

constrast_acf_estimate() {

all_fwhm_files=()

for h in "${hierarchy[@]}";do
	for d in "${duration[@]}";do

		outDir=${outRoot}/${h}_${d}_FWHM
		cd ${outDir}
		
		acf_tmp=$(awk '{a+=$1; b+=$2; c+=$3; n++} END {print a/n, b/n, c/n}' group_ACF_ests.1D) 
		outFile=${outRoot}/${h}_${d}_FWHM.1D
		echo "${acf_tmp}" > ${outRoot}/${h}_${d}_FWHM.1D> ${outFile}
		echo "++ The group average ACF params are: ${acf_tmp}"
		
		all_fwhm_files+=("${outFile}")
		
	done
done

for i in "${!contrasts[@]}";do
	outDir=${outRoot}/ClustSim_${contrasts[$i]}
	mkdir -p ${outDir}
	cd ${outDir}
	
	index_list="${contrast_indices[$i]}"
	for idx in $index_list; do
	echo $idx
      	  files+=("${all_fwhm_files[$idx]}")
    	done
	
	  # Average ACF values using awk
    acf_avg=$(awk '{a+=$1; b+=$2; c+=$3; n++} END {print a/n, b/n, c/n}' ${files[@]})

    # Save and report
    echo $acf_avg > ${outDir}/${contrasts[$i]}_avgFWHM.1D
    echo "++ Saved average ACF for:" ${contrasts[$i]}
   
done

}

cluster_simulation(){

for con in "${contrasts[@]}";do
	outDir=${outRoot}/ClustSim_${con}
	cd ${outDir}
	
	acf_params=$(<${outDir}/${con}_avgFWHM.1D)

	3dClustSim -both -mask ${maskFile} -acf ${acf_params} -prefix ClustSim 
done
}


get_clust_thresh () {
  local ClustSimPath="$1"
  local nn="$2"
  local direct="$3"
  local pVal="$4"
  local alphaPos="$5"

  local clustsim_file="${ClustSimPath}/ClustSim.NN${nn}_${direct}.1D"
  #echo "Looking for p-value: $pVal at column $alphaPos"
  awk -v pval="$pVal" -v pos="$alphaPos" '
    $1 == pval { print $pos }
  ' "$clustsim_file"
}


clusterize_map() {

for con in "${contrasts[@]}";do
	
	conDir=${studyDir}/_AFNI_Analysis/${con}
	clustOutDir=${conDir}/ClusterizeOutput_report
	mkdir -p ${clustOutDir}

	cd ${clustOutDir}
	
	outDir=${outRoot}/ClustSim_${con}

	voxThresh=$(get_clust_thresh ${outDir} ${neighbors} ${testDirect} ${pThresh} 7) # alpha 0.05 is at position 7
	echo "Voxel threshold for p=$pThresh and alpha=0.05: $voxThresh"
	
	textOut=${clustOutDir}/${con}_report.1D


	3dClusterize -inset ${conDir}/Con_${con}+orig \
		-ithr 1 \
		-idat 0 \
		-mask ${maskFile} \
		-NN ${neighbors} \
		-bisided p=${pThresh} \
		-clust_nvox ${voxThresh} \
		-pref_map Con_${con}_Map \
		-pref_dat Con_${con}_EffEst \
		-orient LPI \
		-1Dformat > ${textOut}
		
	# Split Output into Positive and negative 
	3dcalc -a Con_${con}_Map+orig -b Con_${con}_EffEst+orig \
		-expr 'a*step(b)' -prefix Pos_${con}_Map
	3dcalc -a Con_${con}_Map+orig -b Con_${con}_EffEst+orig \
		-expr 'a*step(-b)' -prefix Neg_${con}_Map	
		
	3dcalc -a Con_${con}_EffEst+orig  -expr 'a*step(a)' -prefix Pos_${con}_EffEst
	3dcalc -a Con_${con}_EffEst+orig  -expr '-a*step(-a)' -prefix Neg_${con}_EffEst
	
	3dAFNItoNIFTI -prefix Con_${con}_Map.nii Con_${con}_Map+orig
	3dAFNItoNIFTI -prefix Con_${con}_EffEst.nii Con_${con}_EffEst+orig
	
	3dAFNItoNIFTI -prefix Pos_${con}_Map.nii Pos_${con}_Map+orig
	3dAFNItoNIFTI -prefix Pos_${con}_EffEst.nii Pos_${con}_EffEst+orig
	3dAFNItoNIFTI -prefix Neg_${con}_Map.nii Neg_${con}_Map+orig
	3dAFNItoNIFTI -prefix Neg_${con}_EffEst.nii Neg_${con}_EffEst+orig


done

}


localize_clusters () {

for con in "${contrasts[@]}";do
	
	conDir=${studyDir}/_AFNI_Analysis/${con}
	clustOutDir=${conDir}/ClusterizeOutput_report
	cd ${clustOutDir}
	
	textOut=${clustOutDir}/${con}_report.1D
	
	localOut=${clustOutDir}/${con}_peaklocations.text
	
	whereami -coord_file ${textOut}'[1,2,3]' -tab -lpi -space MNI > ${localOut} || true

done

}


local_extrema () {
for con in "${contrasts[@]}";do
	echo ${con}
	conDir=${studyDir}/_AFNI_Analysis/${con}
	clustOutDir=${conDir}/ClusterizeOutput_report
	cd ${clustOutDir}
	
	binClustFile=Con_${con}_Map+orig

	read min max <<< $(3dBrickStat -min -max ${binClustFile})
	echo "Min: $min"
	echo "Max: $max"
	
	right_thr=$(grep 'Threshold value' ${con}_report.1D | awk -F'thr=' '{print $3}' | awk -F';' '{print $1}'| tr -d ']')

	echo ${right_thr}
	rm -f tstat_tmp+orig* temp_mask+orig*
	
	if [ "$max" -gt 0 ]; then
		for clustID in $(seq 1 $max); do
			echo $clustID
			3dcalc -a ${binClustFile} -expr "equals(a,${clustID})" -prefix temp_mask
			3dcalc -a temp_mask+orig -b ${conDir}/Con_${con}+orig'[1]' -expr "a*b" -prefix tstat_tmp
			max_t=$(3dBrickStat -max tstat_tmp+orig | awk '{print $1}')
			echo "Max: $max_t"
			if (( $(echo "$max_t < $right_thr" | bc -l) )); then
				thresh=$(echo "-1 * $right_thr" | bc -l)
				echo "$thresh"
				3dExtrema -mask_file temp_mask+orig \
					  -prefix Con_${con}_cluster_${clustID}_peaks \
					  -data_thr "$thresh" \
					  -nbest 5 \
					  -sep_dist 8 \
					  -minima -volume "${conDir}/Con_${con}+orig[1]"
			else
				thresh=$(echo $right_thr)
				echo "$thresh"
				3dExtrema -mask_file temp_mask+orig \
					  -prefix Con_${con}_cluster_${clustID}_peaks \
					  -data_thr "$thresh" \
					  -nbest 5 \
					  -sep_dist 8 \
					  -maxima -volume "${conDir}/Con_${con}+orig[1]"
			fi
			
		
			3dmaskdump -mask Con_${con}_cluster_${clustID}_peaks+orig -xyz -noijk Con_${con}_EffEst+orig | awk '{print -$1, -$2, $3, $4}' \
			> Con_${con}_cluster_${clustID}_peaks_coords.txt

			
			rm tstat_tmp+orig* temp_mask+orig*
						
		done
	fi
		
done
}

# Run functions
estimate_smoothness
constrast_acf_estimate
cluster_simulation
clusterize_map
localize_clusters
local_extrema



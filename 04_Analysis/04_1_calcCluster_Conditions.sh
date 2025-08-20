#!/bin/bash
set -e

# Run3dIsc_MH_v02.sh
#
# Created 2025-04-23 by FM.


# Variables
studyDir=/.../NIFTI


subjects=("BHM_437")

hierarchy=("Shot" "Scene")
duration=("4s" "12s" "36s")

maskFile=${studyDir}/bin_GM_mask.nii

outRoot=${studyDir}/_AFNI_Analysis/00_Conditions/ISC3D_Con
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

for h in "${hierarchy[@]}";do
	for d in "${duration[@]}";do

		outDir=${outRoot}/${h}_${d}_FWHM
		cd ${outDir}
		
		acf_tmp=$(awk '{a+=$1; b+=$2; c+=$3; n++} END {print a/n, b/n, c/n}' group_ACF_ests.1D) 
		outFile=${outDir}/${h}_${d}_FWHM.1D
		echo "${acf_tmp}" > ${outDir}/${h}_${d}_FWHM.1D> ${outFile}
		echo "++ The group average ACF params are: ${acf_tmp}"
		
	done
done

}

cluster_simulation(){

for h in "${hierarchy[@]}";do
	for d in "${duration[@]}";do

		outDir=${outRoot}/${h}_${d}_FWHM
		cd ${outDir}
	
		acf_params=$(<${outDir}/${h}_${d}_FWHM.1D)

		3dClustSim -both -mask ${maskFile} -acf ${acf_params} -prefix ClustSim 
	done
done
}


get_clust_thresh () {

  clustsim_file=${outRoot}/${h}_${d}_FWHM/ClustSim.NN${neighbors}_${testDirect}.1D
  # alpha 0.05 is at position 7 in that line

awk -v pval="$pThresh" '
    $1 == pval { print $7 } 
  ' "$clustsim_file"
  
}


clusterize_map() {

clustOutDir=${outRoot}/ClusterizeOutput
mkdir -p ${clustOutDir}

cd ${clustOutDir}	

for h in "${hierarchy[@]}";do
	for d in "${duration[@]}";do

		voxThresh=$(get_clust_thresh)
		echo "Voxel threshold for p=$pThresh and alpha=0.05: $voxThresh"


		3dClusterize -inset ${outRoot}/${h}_${d}+orig \
			-ithr 1 \
			-idat 0 \
			-mask ${maskFile} \
			-NN ${neighbors} \
			-bisided p=${pThresh} \
			-clust_nvox ${voxThresh} \
			-pref_map ${h}_${d}_Map \
			-pref_dat ${h}_${d}_EffEst
		
		3dAFNItoNIFTI -prefix ${h}_${d}_Map.nii ${h}_${d}_Map+orig
		3dAFNItoNIFTI -prefix ${h}_${d}_EffEst.nii ${h}_${d}_EffEst+orig

	done
done

}


# Run functions
estimate_smoothness
constrast_acf_estimate
cluster_simulation
clusterize_map



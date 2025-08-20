#!/bin/bash

# Create all necessary contrast maps to run AFNI 3dISC analysis according to example 5/ example 1

# for that idea we calcualte the difference between the pairwise correlation maps while additionally z-scaling them/or I directly use the z-scaled ones

# in total we have 13 comparisions we need to cacualte
# 1.Shot vs Scene
# 2. 04s vs 12s
# 3. 04s vs 36s
# 4. 12s vs 36s
# 5. Shot04s vs Scene04s
# 6. Shot12s vs Scene12s
# 7. Shot36s vs Scene36s
# 8. Shot04s vs Shot12s
# 9. Shot04s vs Shot36s
# 10. Shot12s vs Shot36s
# 11. Scene04s vs Scene12s
# 12. Scene04s vs Scene36s
# 13. Scene12s vs Scene36s

# Different effect of Scene vs. Shot in different levels of Duration
# 1.(Shot04s vs Scene04s) vs (Shot12s vs Scene12s)
# 2.(Shot04s vs Scene04s) vs (Shot36s vs Scene36s)
# 3.(Shot12s vs Scene12s) vs (Shot36s vs Scene36s)

# We just need the contrasts in one direction bcause we will check the contrast bi-directionally to get the effects in both directions

# Variables
studyDir=/.../NIFTI

subjects=("BHM_347")


# 1.Shot vs Scene
Shot_Scene() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"01_Shot_Scene"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-c ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-d ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-e ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-f ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr '((1/3)*atanh(d) + (1/3)*atanh(e) + (1/3)*atanh(f)) - ((1/3)*atanh(a) + (1/3)*atanh(b) + (1/3)*atanh(c))' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 2. 04s vs 12s
d04s_d12s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"02_04s_12s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-c ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-d ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-expr '((1/2)*atanh(b) + (1/2)*atanh(d)) - ((1/2)*atanh(a) + (1/2)*atanh(c))' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 3. 04s vs 36s
d04s_d36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"03_04s_36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-c ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-d ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr '((1/2)*atanh(b) + (1/2)*atanh(d)) - ((1/2)*atanh(a) + (1/2)*atanh(c))' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 4. 12s vs 36s
d12s_d36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"04_12s_36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-c ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-d ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr '((1/2)*atanh(b) + (1/2)*atanh(d)) - ((1/2)*atanh(a) + (1/2)*atanh(c))' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 5. Shot04s vs Scene04s
Shot04s_Scene04s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"05_Shot04_Scene04"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 6. Shot12s vs Scene12s
Shot12s_Scene12s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"06_Shot12s_Scene12s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 7. Shot36s vs Scene36s
Shot36s_Scene36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"07_Shot36s_Scene36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 8. Shot04s vs Shot12s
Shot04s_Shot12s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"08_Shot04s_Shot12s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 9. Shot04s vs Shot36s
Shot04s_Shot36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"09_Shot04_Shot36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 10. Shot12s vs Shot36s
Shot12s_Shot36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"10_Shot12s_Shot36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_12s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Shot_36s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 11. Scene04s vs Scene12s
Scene04s_Scene12s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"11_Scene04s_Scene12s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 12. Scene04s vs Scene36s
Scene04s_Scene36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"12_Scene04s_Scene36s"
	mkdir -p ${outDir}
	cd ${outDir}
	echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_4s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 13. Scene12s vs Scene36s
Scene12s_Scene36s() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"13_Scene12s_Scene36s"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_12s_"${subj1}_${subj2}.nii \
				-b ${subjDir}/"ISC_Pair_maps_M13"/"ISC_Pair_corr_Scene_36s_"${subj1}_${subj2}.nii \
				-expr 'atanh(b) - atanh(a)' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 1.(Shot04s vs Scene04s) vs (Shot12s vs Scene12s) bzw. (Scene12s - Shot12s) - (Scene04s - Shot04s)
INT01_12sSceneShot_04sSceneShot() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"INT01_12sSceneShot_04sSceneShot"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"06_Shot12s_Scene12s"/${subj1}_${subj2}+orig \
				-b ${subjDir}/"05_Shot04_Scene04"/${subj1}_${subj2}+orig \
				-expr 'a-b' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 2.(Shot04s vs Scene04s) vs (Shot36s vs Scene36s) bzw. (Scene36s - Shot36s) - (Scene04s - Shot04s)
INT02_36sSceneShot_04sSceneShot() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"INT02_36sSceneShot_04sSceneShot"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"07_Shot36s_Scene36s"/${subj1}_${subj2}+orig \
				-b ${subjDir}/"05_Shot04_Scene04"/${subj1}_${subj2}+orig \
				-expr 'a-b' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

# 3.(Shot12s vs Scene12s) vs (Shot36s vs Scene36s) bzw. (Scene36s - Shot36s) - (Scene12s - Shot12s)
INT03_36sSceneShot_12sSceneShot() {
for subj1 in "${subjects[@]}";do
	
	subjDir=${studyDir}/${subj1}
	outDir=${subjDir}/"INT03_36sSceneShot_12sSceneShot"
	mkdir -p ${outDir}
	cd ${outDir}
	 echo "calcualting contrasts for ${subj1}"

	for subj2 in "${subjects[@]}";do
	
		if [ "$subj1" != "$subj2" ]; then
        	    # Only runs if subj1 and subj2 are different
       	    
	            	# your command here, e.g., some_command "$subj1" "$subj2"
			3dcalc -a ${subjDir}/"07_Shot36s_Scene36s"/${subj1}_${subj2}+orig \
				-b ${subjDir}/"06_Shot12s_Scene12s"/${subj1}_${subj2}+orig \
				-expr 'a-b' \
				-prefix ${subj1}_${subj2}
		fi
	done
done
}

############################################

# RUN COMMANDS

# Main effects Hierarchy
time Shot_Scene

# Main effect Duration
d04s_d12s
d04s_d36s
d12s_d36s

# Effect of Hierarchy in Levels of Duration
Shot04s_Scene04s
Shot12s_Scene12s
Shot36s_Scene36s

# Effect of Duration in Levels of Hierarchy
Shot04s_Shot12s
Shot04s_Shot36s
Shot12s_Shot36s
#
Scene04s_Scene12s
Scene04s_Scene36s
Scene12s_Scene36s

# Different effect of Scene vs. Shot in different levels of Duration
time INT01_12sSceneShot_04sSceneShot 
INT02_36sSceneShot_04sSceneShot
INT03_36sSceneShot_12sSceneShot

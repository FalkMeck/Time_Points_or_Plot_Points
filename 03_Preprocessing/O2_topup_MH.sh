# FSL Setup
FSLDIR=/.../fsl
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh

# My custom values for FSL environment variables
export FSLOUTPUTTYPE=NIFTI_GZ

# declare -a subjects=("ABN_653" "NAN_439" "MIN_315" "GGH_758" "DMH_110") # EXAMPLE
declare -a subjects=()


ANALYSISDIR=/...

ACQP=$ANALYSISDIR/acqparams.txt

NiftiDIR=/.../NIFTI

for subj in "${subjects[@]}"
do
	echo $subj
	SUBJECTDIR="$NiftiDIR"/$subj/topup
	
	cd "$SUBJECTDIR"
	
	INFILE=4Dfortopup

	topup --imain="$INFILE" --datain="$ACQP" --config=b02b0_1.cnf --out=topup_results \
		--fout=field_names --iout=unwarped_images
		
	for file in epis*; do
    		# Extract the filename without the directory and extension
    		found_file=$(basename "${file%.*}")
    		outName="u_${found_file}"
    		echo "$outName"
   		# Run the applytopup command
    		applytopup --imain="$found_file" --topup=topup_results --datain="$ACQP" --inindex=1 \
             		  --out="$outName" --method=jac
	done
		
	echo "DONE"
done		

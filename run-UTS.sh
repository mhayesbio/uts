#!/bin/bash
# Goal: Condense commands to run MMQuant followed by UTS script
# MAKE PYTHON EXECUTABLE (BEFORE USING BASH job.sh): chmod u+x job.sh
# Optional: run DGE analysis using R DESeq2

# Path to mmquant
PATH_MMQUANT=$(which mmquant)

# Path to samtools
PATH_SAMTOOLS=$(which samtools)

# Path to python v3 or greater
PATH_PYTHON=$(which python)

# Path to Rscript (optional, not required for UTS correction)
PATH_RSCRIPT=$(which Rscript)

#####################################################

# Ensure that mmquant file exists
if [[ ! -f $PATH_MMQUANT ]] || [[ $(basename $PATH_MMQUANT) != mmquant ]]; then
	printf "\nmmquant does not exist in path \"$PATH_MMQUANT\"."
	printf "\nCheck path to mmquant file is defined by PATH_MMQUANT.\nAborting.\n\n"
    	exit
fi

# Check samtools file exists
if [[ ! -f $PATH_SAMTOOLS ]] || [[ $(basename $PATH_SAMTOOLS) != samtools ]]; then
	printf "\nsamtools does not exist in path \"$PATH_SAMTOOLS\"."
	printf "\nCheck path to samtools file is defined by PATH_SAMTOOLS.\nAborting.\n\n"
	exit
fi

# Check python file exists
if [[ ! -f $PATH_PYTHON ]]; then
	printf "\npython does not exist in path \"$PATH_PYTHON\"."
	printf "\nCheck path to python file is defined by PATH_PYTHON.\nAborting.\n\n"
	exit
fi
if [[ $(basename $PATH_PYTHON) != python ]] && [[ $(basename $PATH_PYTHON | rev | cut -c 2- | rev) != python ]]; then
	printf "\npython does not exist in path \"$PATH_PYTHON\"."
	printf "\nCheck path to python file is defined by PATH_PYTHON.\nAborting.\n\n"
	echo "TRUE"
	exit	
fi

# Check python version is 3 or greater
if [[ $($PATH_PYTHON -V | cut -d ' '  -f2 | cut -c 1) -lt 3 ]]; then
	printf "\nPython Version 3 or greater is required. Current version is $($PATH_PYTHON -V | cut -d ' ' -f2).\nAborting.\n\n"
	exit
fi

# Get parameters
ERCC='FALSE'
while getopts a:s:b:o:e:m:d:f:-: flag; do
	case "${flag}" in
		a) GTF=${OPTARG};;
		s) STRAND=${OPTARG};;
		b) BAMS=${OPTARG};;
		o) OUT=${OPTARG};;
		e) EFF=${OPTARG};;
		m) DGE_MIX=${OPTARG};;
		d) DILUTION=${OPTARG};;
		f) STM=${OPTARG};;
		-)
			case "${OPTARG}" in
				ercc) ERCC='TRUE';;
				mpc) MPC='TRUE';;
				dge) DGE='TRUE';;
				disable) DISABLE='TRUE';;
			esac	
	esac
done

# CHeck required  parameters are defined
if [[ -z $GTF ]] || [[ -z $STRAND ]] || [[ -z $BAMS ]] || [[ -z $OUT ]]; then
	printf "\nRequired parameters\n"
	printf "\t-a     <GTF file>\n"
	printf "\t-s     <R | F> strandedness\n"
	printf "\t-b     <BAM file | directory to BAM files>\n"
	printf "\t-o     <output directory>\n"
	printf "\nOptional paramters\n"
	printf "\t-e     <integer> mean fragment length (default = 200)\n"
	printf "\t--mpc  Molecules Per Cell quantificaiton using ERCC Spike-In RNAs. Requires -m and -d\n"
	printf "\t-m     <1 | 2> ERCC Mix 1 or Mix 2. Used with --mpc and --ercc\n"
	printf "\t-d     <integer> ERCC dilution factor per 1e5 cells. Used with --mpc\n"
	printf "\t--dge  Differential Gene Expression  analysis. Requires -f\n"
	printf "\t-f     <Sample Treatment Matrix file>. Used with --dge\n"
	printf "\t--ercc ERCC normalized Differential Gene Expression analysis. Used with --dge. Can use -m to specify which ERCC Spike-In RNA mix to normalize by. Default is to normalize by ERCC Spike-In RNA group B genes which contain the same concentrations in both Mix 1 and Mix 2.\n"
	printf "\t--disable Disables mmquant and UTS correction. Useful for running other analyses if UTS output has already been generated.\n\n"
	exit
fi

# Checks
if [[ $DISABLE != TRUE ]]; then
	# Check GTF file
	if [[ ! -f $GTF ]] || [[ $(basename $GTF | rev | cut -c -3 | rev) != gtf ]]; then
		printf "\nError: incorrect GTF file input.\nAborting.\n\n"
		exit
	fi
	
	# Check STRAND inuput is correct
	if [[ $STRAND != R ]] && [[ $STRAND != F ]]; then
		printf "\nError: incorrect library strandedness input.\nAborting.\n\n"
		exit
	fi

	# Check output directory exists
	if [[ ! -d ${OUT} ]]; then
		printf "\nError: output directory \"$OUT\" does not exist.\nAborting.\n\n"
		exit
	fi

fi

# If MPC enabled, check -d and -m are specified
if [[ $MPC == TRUE ]]; then
	# Check dilution and mix input specified
	if [[ -z $DGE_MIX ]] && [[ -z $DILUTION ]]; then
		printf "\nError: --mpc requires -m and -d options.\nAborting.\n\n"
		exit
	fi
	
	# Check Mix is specified
	if [[ $DGE_MIX -ne 1 ]] && [[ $DGE_MIX -ne 2 ]]; then
		printf "\nError: -m must equal 1 or 2.\nAborting.\n\n"
		exit
	fi
	
	# Check dilution factor is a numerical value
	if [[ "$DILUTION" =~ ^[+-]?[0-9]+$ ]] || [[ $DILUTION =~ ^[+-]?[0-9]*\.[0-9]+$ ]]; then
		:
	else
		printf "\nError: -d must be a numerical value. Non-integers are allowed.\n\n"
		exit
	fi
fi

# If DGE enabled, check -f specified and if -m is specified
if [[ $DGE == TRUE ]]; then
	if [[ -z $STM ]]; then
		printf "\nError: --dge requires -f option.\nAborting.\n\n"
		exit
	fi

	if [[ ! -f $STM ]]; then
		printf "\nError: -f option must be a file.\nAborting.\n\n"
		exit
	fi

	# Check DGE_MIX input
	if [[ ! -z $DGE_MIX ]]; then
	       if [[ $DGE_MIX -ne 1 ]] && [[ $DGE_MIX -ne 2 ]] && [[ $DGE_MIX != group_b ]]; then
		       printf "\nError: incorrect -m input..\nAborting\n\n"
		       exit
	       fi
	fi

	# If ERCC Spike-In RNAs Mix 1 or Mix 2 not specifeid, set normalization to group B genes
	if [[ $DGE_MIX -ne 1 ]] && [[ $DGE_MIX -ne 2 ]]; then
		DGE_MIX="group_b"
	fi
fi

# Check BAM files exist
if [[ $(find $BAMS -maxdepth 1 -type f -name "*.bam" -not -name "*Aligned.toTranscriptome.out.bam" -exec ls {} + | wc -l) -eq 0 ]]; then
	printf "\nError: no .bam file(s) detected in \"$PATH_BAM\".Aborting.\n\n"
	exit
fi

# Get list of BAM files
BAMSLIST=$( find $BAMS -maxdepth 1 -type f -name "*.bam" -not -name "*.Aligned.toTranscriptome.out.bam" -exec ls {} + )
# Get array of BAM files
BAMSARRAY=($(ls -d $BAMSLIST))

# Get NCUT value for unique file name 
if [[ $(echo "$BAMSLIST" | grep .Aligned.sortedByCoord.out.bam | wc -l) -gt 0 ]]; then
	NCUT=31
else
	if [[ $(echo "$BAMSLIST" | grep .Aligned.out.bam | wc -l) -gt 0 ]]; then
		NCUT=17
	else
		NCUT=5
	fi
fi

# Create output directory
if [[ ! -d ${OUT}/UTS_output ]]; then
	mkdir ${OUT}/UTS_output
fi

STARTDATE=$(date)
SECONDS=0

# Generate mmquant files
# COMMAND: /mmquant -a annotation.gtf -r sample1.bam -o <output_file_name> -s <strandedness (R/F)> -f BAM
if [[ $DISABLE != TRUE ]]; then
	printf "\nRunning mmquant...\n"
	for BAM in $BAMSLIST; do
		SAVE_SECONDS=$SECONDS
		SECONDS=0
		BAMNAME=$( echo $BAM | sed 's!.*/!!' | rev | cut -c $NCUT- | rev | sed 's!.*/!!' )
		$PATH_MMQUANT -a $GTF -r $BAM -o ${OUT}/UTS_output/${BAMNAME}_temp -s $STRAND -f BAM
		DURATION=$SECONDS
		SECONDS=$(($SAVE_SECONDS+$DURATION))
		printf "$BAMNAME done. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n\n"
	done
fi


# Calculate sequencing depth
if [[ $DISABLE != TRUE ]] || [[ $DGE == TRUE ]] || [[ $MPC == TRUE ]]; then
	printf "\nCalculating sequencing depth...\n"
	for BAM in $BAMSLIST; do
		SAVE_SECONDS=$SECONDS
		SECONDS=0
		BAMNAME=$( echo $BAM | sed 's!.*/!!' | rev | cut -c $NCUT- | rev | sed 's!.*/!!' )
		DEPTH+=("$($PATH_SAMTOOLS view -c -q 255 $BAM)")
		DURATION=$SECONDS
		SECONDS=$(($SAVE_SECONDS+$DURATION))
		printf "$BAMNAME done. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n"
	done
	printf "\n"
fi

# Run UTS script
# COMMAND: Effective_Length_UTS.py <GTF-file> <mmquant-file> <output-file-name> <samtools scaling factor> <optional: mean fragment length>
# Samtools scaling factor via: samtools view -c -q 255 file.bam
# Ref: https://stackoverflow.com/questions/9332802/how-to-write-a-bash-script-that-takes-optional-input-arguments
if [[ $DISABLE != TRUE ]]; then
	printf "\nRunning UTS corrections...\n"
	for ((i=0;i<=$(($(ls $BAMSLIST | wc -l)-1));i++)); do
		SAVE_SECONDS=$SECONDS
		SECONDS=0
		BAMNAME=$(echo "${BAMSARRAY[$i]}" | sed 's!.*/!!' | rev | cut -c $NCUT- | rev | sed 's!.*/!!')
		if [[ -f ${OUT}/UTS_output/${BAMNAME}.genes.results ]]; then
			rm ${OUT}/UTS_output/${BAMNAME}.genes.results
		fi
		$PATH_PYTHON $(dirname $0)/bin/run_UTS.py $GTF ${OUT}/UTS_output/${BAMNAME}_temp ${OUT}/UTS_output/${BAMNAME}.genes.results "${DEPTH[$i]}" $EFF &>/dev/null
		rm ${OUT}/UTS_output/${BAMNAME}_temp
		DURATION=$SECONDS
		SECONDS=$((SAVE_SECONDS+$DURATION))
		printf "$BAMNAME done. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n"
	done
	printf "\n"
fi

# Molecules Per Cell quantitation
if [[ $MPC == TRUE ]]; then
	printf "\nRunning Molecules Per Cell quantitation...\n"
	for ((i=0;i<=$(($(ls $BAMSLIST | wc -l)-1));i++)); do
		SAVE_SECONDS=$SECONDS
		SECONDS=0
		BAMNAME=$(echo "${BAMSARRAY[$i]}" | sed 's!.*/!!' | rev | cut -c $NCUT- | rev | sed 's!.*/!!')
		$PATH_RSCRIPT $(dirname $0)/bin/ERCCnorm_UTS.R ${OUT}/UTS_output/${BAMNAME}.genes.results "${DEPTH[$i]}" $DILUTION $DGE_MIX
		DURATION=$SECONDS
		SECONDS=$((SAVE_SECONDS+$DURATION))
		printf "$BAMNAME done. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n"

	done
fi

# Differential Gene Expression Analysis
if [[ $DGE == TRUE ]]; then
	if [[ $ERCC == TRUE ]]; then
		ERCC_DGE="ENABLED"
		printf "\nRunning Differential Gene Expression analysis with ERCC Spike-In RNA normalization and conventional normalization\n"
	else
		ERCC_DGE="DISABLED"
		printf "\nRunning Differential Gene Expression analysis with conventional normalization.\n"
	fi
	SAVE_SECONDS=$SECONDS
	SECONDS=0
	BAMNAMES=$( ls $BAMSLIST | sed 's!.*/!!' | rev | cut -c $NCUT- | rev | sed 's!.*/!!')
	$PATH_RSCRIPT $(dirname $0)/bin/DGE_UTS.R ${OUT}/UTS_output "$(printf "%s " "${DEPTH[@]}" | rev | cut -c 2- | rev)" $ERCC_DGE $DGE_MIX $STM "$(printf "%s " $BAMNAMES | rev | cut -c 2- | rev)"
	DURATION=$SECONDS
	SECONDS=$((SAVE_SECONDS+$DURATION))
	printf "Differential Gene Expression analysis done. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n"
fi

ENDDATE=$(date)
DURATION=$SECONDS

printf "\nDone. Elapsed time: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.\n"

# Annotate UTS_output.Log file
if [[ $DISABLE == TRUE ]]; then
	DISABLE_OUTPUT="mmquant and UTS correction were disabled."
else
	DISABLE_OUTPUT="mmquant and UTS correction were performed."
fi

if [[ $MPC == TRUE ]]; then
	MPC_OUTPUT="Molecules Per Cell quantitation was performed with dilution factor of \"$DILUTION\" per 1e5 cells and standardized curve using ERCC Spike-In RNAs from Mix \"$DGE_MIX\""
else
	MPC_OUTPUT="Molecules Per Cell quantitation was disabled"
fi

if  [[ $DGE == TRUE ]]; then
	if  [[ $ERCC == TRUE ]]; then
		if [[ $DGE_MIX -ne 1 ]] && [[ $DGE_MIX -ne 2 ]]; then
			DGE_OUTPUT="Differential Gene Expression analysis was performed usig ERCC Spike-In RNA group B genes and conventional normalization."
		else
			DGE_OUTPUT="Differential Gene Expression analysis was performed using all ERCC Spike-In RNA genes."
		fi
	else
		DGE_OUTPUT="Differential Gene Expression analysis was performed using conventional normalization."
	fi
else
	DGE_OUTPUT="Differential Gene Expression analysis was disabled"
fi

# Annotate .Log file
cat <<EOT >> ${OUT}/temp
START DATE: $STARTDATE
END DATE: $ENDDATE
Duration of run: $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds.

GTF file used: $(basename $GTF)
Strandedness: $STRAND
$DISABLE_OUTPUT
$MPC_OUTPUT
$DGE_OUTPUT

BAM files analyzed: 
$(find $BAMS -maxdepth 1 -type f -name "*bam" -not -name "*Aligned.toTranscriptome.out.bam" -printf '%f\n')
EOT

if [[ -f ${OUT}/UTS_output.Log ]]; then
	cat << EOT>> ${OUT}/temp

####################################################################

EOT
	cat ${OUT}/UTS_output.Log >> ${OUT}/temp
	rm ${OUT}/UTS_output.Log
	mv ${OUT}/temp ${OUT}/UTS_output.Log
else
	mv ${OUT}/temp ${OUT}/UTS_output.Log
fi



#!/bin/bash
#$ -cwd 


##########################################################################
# Allysia Mak                                                            #
# mutect2.sh runs SplitNCigarReads for read preprocessing before running #
# Mutect2 on all bams in the 'parentdir'                                 #
#                                                                        #
# The number of bams run in parallel is specified by 'numProcesses'      #
# Input to stdin: $1 = number of processes $2 = 'pon' mode or null for   #
#     matched normal mode                                                # 
# Change lines 21-27 to fit your file structures                         # 
#                                                                        #
# !! Make sure there is only 1 '*.bam' file per sample, if more than one #
#    then change first for loop to find -name '*your.bam'                #
##########################################################################

# Change parentdir to directory running on (normal or tumor bam directory)

#parentdir="/scratch/jakutagawa/RNA-seq/GTEx_bams/STAR_aligned/subset-SRA"
parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/tumor"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/normal"
hs37="/pod/pstore/groups/brookslab/reference_indices/hs37/hs37d5.fa"
hg19="/pod/pstore/groups/brookslab/reference_indices/hg19/hg19.fa"
outDir="/scratch/amak/varCalls/Mutect/pon-tumor-vcfs"
pon="/pod/pstore/groups/brookslab/amak/var-prediction-testing/gtex-vcfs/GTEX_PON.vcf"
numProcesses=$1
matchedOrPON=$2
if [ $((numProcesses)) -eq 0 ]; then
    echo "Max number of processes not supplied. Defaulting to 2"
    numProcesses=2
fi


function maxProcesses {
    # Waits until there are less than 'numProcesses' jobs running before starting a new job

    while [ $(jobs | wc -l) -ge $numProcesses ]; do
	#echo 'waiting'
	sleep 5

    done
}

##### Comment this for loop out if split bams are already made by SplitNCigarReads for all bams
for bam in $(find $parentdir -mindepth 1 -name '*.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    uid=$(echo $1)
    unset IFS
    bamdir=$(dirname $bam)

    if [[ -z "$(find $parentdir -name $uid'.split.bam')" ]]; then 
	echo "SplitNCigarReads running on " $uid

	maxProcesses; nice time gatk SplitNCigarReads -fixNDN -R $hs37 -I $bam -O $bamdir'/'$uid'.split.bam' &
	echo $(find $parentdir -name $uid'.split.bam')
    else
	echo $uid ' already done' 
    fi
    
done
wait

for bam in $(find $parentdir -mindepth 1 -name '*.split.bam'); do  #'*split_picard_edited.bam'); do #'*hs37d5.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    uid=$(echo $1)
    unset IFS
    bamdir=$(dirname $bam)
    if [ ! -d "$outDir/$uid" ]; then
        mkdir $outDir/$uid

	echo $bam
	normal_ID_line=$(samtools view -H $bam | grep '@RG')


	IFS=$':\t'
	set $normal_ID_line
	## !!! May need to change depending on bam header !! ##
	#GTEX bams, normal ID is at $9
#	normal_ID=$(echo $9)
	#PCAWG bams, normal ID is at $6
	normal_ID=$(echo $6)
	unset IFS
    
	echo 'normal_ID ' $normal_ID

	#    maxProcesses; nice time gatk Mutect2 -R $hs37 -I $bamdir'/'$uid'.split.bam' -tumor $normal_ID -O $outDir'/'$uid'/'$uid'.vcf' &
	## if running in no-match (tumor-only) mode, take out -pon argument
	if [[ $matchedOrPON = "pon" ]]; then
	    maxProcesses; nice time gatk Mutect2 -R $hs37 -I $bam -tumor $normal_ID -pon $pon -O $outDir'/'$uid'/'$uid'.vcf' &  
	else
	    maxProcesses; nice time gatk Mutect2 -R $hs37 -I $bam -tumor $normal_ID -O $outDir'/'$uid'/'$uid'.vcf' &
	fi
    fi 
#    rm $tumordir'/'$uid'.split.bam'                                                                                                          

    echo $uid " complete"


done
wait

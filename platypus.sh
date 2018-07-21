#!/bin/bash
#$ -cwd 


##########################################################################
# Allysia Mak                                                            #
# platypus.sh runs Opposum for read preprocessing before running         #
# Platypus on all bams in the 'parentdir'                                #
#                                                                        #
# The number of bams run in parallel is specified by 'numJobs'           #
# Input to stdin: $1 = number of processes                               #                                                                    
# Change lines 18-24 to fit your file structures                         #                                                                    
#                                                                        #                                                                    
# !! Make sure there is only 1 '*.bam' file per sample, if more than one #                                                                    
#    then change first for loop to find -name '*your.bam'                # 
##########################################################################

parentdir="/scratch/jakutagawa/RNA-seq/GTEx_bams/STAR_aligned/subset-SRA"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/tumor"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/normal"
hs37="/pod/pstore/groups/brookslab/reference_indices/hs37/hs37d5.fa"
hg19="/pod/pstore/groups/brookslab/reference_indices/hg19/hg19.fa"
outDir="/scratch/amak/varCalls/Platypus/gtex_normals"
opossum="/pod/pstore/groups/brookslab/amak/packages/Opossum/Opossum.py"
platypus="/pod/pstore/groups/brookslab/amak/packages/Platypus/bin/Platypus.py"
#Specify here the number of concurrent jobs
maxProcesses=$1
if [ $(($maxProcesses)) -eq 0 ]; then
    echo "Max number of processes not supplied. Defaulting to 2"
    maxProcesses=2
fi


function maxJobs {
    # Waits until there are less than 'numJobs' jobs running before starting a new job
    while [ $(jobs | wc -l) -ge $maxProcesses ]; do
#	echo 'waiting'
	sleep 5
    done
}

# Comment out if .opossum bams have already been made
for bam in $(find $parentdir -mindepth 1 -name '*.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    uid=$(echo $1)
    unset IFS
    bamdir=$(dirname $bam)

#    mkdir $outDir/$uid
    if [[ -z "$(find $parentdir -name $uid'.opossum.bam')" ]]; then
	echo "Opossum running on " $uid
	if [[ $uid == *"tophat"* ]]; then
	    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$bamdir'/'$uid'.opossum.bam' &
	else
	    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$bamdir'/'$uid'.opossum.bam' &
	fi
    else
	echo $uid ' already done'
    fi
done
wait


for bam in $(find $parentdir -mindepth 1 -name '*opossum.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    uid=$(echo $1)
    unset IFS
    bamdir=$(dirname $bam)
    if [ ! -d "$outDir/$uid" ]; then
	mkdir $outDir/$uid
	echo "Platypus running on " $uid

	maxJobs; nice time python $platypus callVariants --bamFiles $bam --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $outDir'/'$uid'/'$uid'.vcf' &



	echo $uid " complete"
    fi
done
wait

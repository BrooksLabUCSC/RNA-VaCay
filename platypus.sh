#!/bin/bash
#$ -cwd


##########################################################################
# Allysia Mak and Jon Akutagawa                                          #
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

synthetic_dir="/scratch/jakutagawa/icgc/bams/synthetic/spiked"
normal_dir="/scratch/jakutagawa/icgc/bams/normal"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/tumor"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/normal"
hs37="/private/groups/brookslab/reference_indices/hs37/hs37d5.fa"
hg19="/private/groups/brookslab/reference_indices/hg19/hg19.fa"
outDir="/scratch/jakutagawa/icgc/var_calls/platypus"
opossum="/private/groups/brookslab/jakutagawa/variant_calling/packages/Opossum/Opossum.py"
platypus="/private/groups/brookslab/jakutagawa/variant_calling/packages/Platypus/bin/Platypus.py"

#Specify here the number of concurrent jobs
maxProcesses=$1
if [ $(($maxProcesses)) -eq 0 ]; then
    echo "Max number of processes not supplied. Defaulting to 2"
    maxProcesses=2
fi


#function maxJobs {
    # Waits until there are less than 'numJobs' jobs running before starting a new job
#    while [ $(jobs | wc -l) -ge $maxProcesses ]; do
#	echo 'waiting'
#	sleep 5
#    done
#}

# Comment out if .opossum bams have already been made
for bam in $(find $synthetic_dir -mindepth 1 -name '*.sorted.bam'); do
    base_filename=$(basename $bam)
    echo $base_filename
    matching_normal=$normal_dir/"${base_filename%.with_variants.sorted.bam}".bam
    echo $matching_normal
    IFS='.'
    set $base_filename
    uid=$(echo $2)
    echo $uid
    unset IFS
    bamdir=$(dirname $bam)
    #matching_normal=$(find $normal_dir -mindepth 1 -name '*'$uid'.STAR.v1.bam')

    #if [ ! -d $outDir/$uid ];
    #then
    #    mkdir $outDir/$uid
    #    echo "making folder $outDir/$uid"
    #fi

    echo $bam
    opossum_spiked_file="${bam%.*}".opossum.bam
    opossum_normal_file="${matching_normal%.*}".opossum.bam
    #echo $opossum_spiked_file

    if [ ! -f $opossum_spiked_file ]; then
    echo "Opossum running on" $bam
    #if [[ $uid == *"tophat"* ]]; then
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$bamdir'/'$uid'.opossum.bam' &
        #nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$opossum_spiked_file &
    #else
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$bamdir'/'$uid'.opossum.bam' &
    #nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$opossum_spiked_file
    #fi
    #echo $bamdir'/'$uid'.opossum.bam'
    else
    echo 'opossum already completed on spiked '$uid
    echo $opossum_spiked_file
    fi

    if [ ! -f $opossum_normal_file ]; then
    echo "Opossum running on" $matching_normal
    #if [[ $uid == *"tophat"* ]]; then
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$bamdir'/'$uid'.opossum.bam' &
        #nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$opossum_spiked_file &
    #else
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$bamdir'/'$uid'.opossum.bam' &
    #nice time python $opossum --BamFile=$matching_normal --SoftClipsExist=True --OutFile=$opossum_normal_file
    #fi
    #echo $bamdir'/'$uid'.opossum.bam'
    else
    echo 'opossum already completed on normal' $uid
    echo $matching_normal
    echo $opossum_normal_file
    fi

    platypus_spiked_file=$outDir'/'$uid'.spiked.vcf'
    platypus_normal_file=$outDir'/'$uid'.normal.vcf'

    if [ ! -f $platypus_spiked_file ]; then
        echo "Platypus running on" $opossum_spiked_file
        nice time python $platypus callVariants --bamFiles $opossum_spiked_file --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $platypus_spiked_file
        echo $platypus_spiked_file" completed"
    else
        echo $uid 'spiked file already completed'
    fi

    if [ ! -f $platypus_normal_file ]; then
        echo "Platypus running on" $opossum_normal_file
        nice time python $platypus callVariants --bamFiles $opossum_normal_file --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $platypus_normal_file
        echo $platypus_normal_file" completed"
    else
        echo $uid 'normal file already completed'
    fi


done



#for bam in $(find $synthetic_dir -mindepth 1 -name '*opossum.bam'); do
#    uid=$(basename $bam)
#    IFS='.'
#    set $uid
#    uid=$(echo $2)
#    unset IFS
#    bamdir=$(dirname $bam)
#    #if [ ! -d "$outDir/$uid" ]; then
#    #mkdir $outDir/$uid
#    echo "Platypus running on" $bam
    #echo $bam
#    matching_normal=$(find $normal_dir -mindepth 1 -name '*'$uid'.STAR.v1.bam')
#    echo 'matching normal is' $matching_normal

    #maxJobs; nice time python $platypus callVariants --bamFiles $bam --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $outDir'/'$uid'/'$uid'.vcf' &
#    nice time python $platypus callVariants --bamFiles $matching_normal --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $outDir'/'$uid'.normal.vcf'
#    nice time python $platypus callVariants --bamFiles $bam --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $outDir'/'$uid'.spiked.vcf'

#    echo $outDir'/'$uid'/'$uid'.vcf'



#    echo $uid "complete"
    #fi
#done

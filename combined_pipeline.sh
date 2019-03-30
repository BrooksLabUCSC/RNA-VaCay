#!/bin/bash
#$ -cwd


##########################################################################
# Jon Akutagawa                                                          #
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
pon="/scratch/jakutagawa/icgc/bams/gtex/pon/pon.vcf.gz"
opossum="/private/groups/brookslab/jakutagawa/variant_calling/packages/Opossum/Opossum.py"
platypus="/private/groups/brookslab/jakutagawa/variant_calling/packages/Platypus/bin/Platypus.py"
gatk="/private/groups/brookslab/jakutagawa/tools/gatk-4.0.11.0/gatk"
vardict="/private/groups/brookslab/jakutagawa/variant_calling/packages/VarDictJava/build/install/VarDict/bin/VarDict"

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
    uid=$(basename $bam)
    echo $uid
    IFS='.'
    set $uid
    uid=$(echo $2)
    echo $uid
    unset IFS
    bamdir=$(dirname $bam)
    matching_normal=$(find $normal_dir -mindepth 1 -name '*'$uid'.STAR.v1.bam')

    #if [ ! -d $outDir/$uid ];
    #then
    #    mkdir $outDir/$uid
    #    echo "making folder $outDir/$uid"
    #fi

    echo $bam

################# platypus section
    opossum_spiked_file="${bam%.*}".opossum.bam
    opossum_normal_file="${matching_normal%.*}".opossum.bam
    #echo $opossum_spiked_file

    if [ ! -f $opossum_spiked_file ]; then
    echo "Opossum running on" $opossum_spiked_file
    #if [[ $uid == *"tophat"* ]]; then
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$bamdir'/'$uid'.opossum.bam' &
        #nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$opossum_spiked_file &
    #else
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$bamdir'/'$uid'.opossum.bam' &
    nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$opossum_spiked_file
    #fi
    #echo $bamdir'/'$uid'.opossum.bam'
    else
    echo 'spiked '$uid 'already done'
    fi

    if [ ! -f $opossum_normal_file ]; then
    echo "Opossum running on" $opossum_normal_file
    #if [[ $uid == *"tophat"* ]]; then
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$bamdir'/'$uid'.opossum.bam' &
        #nice time python $opossum --BamFile=$bam --SoftClipsExist=False --OutFile=$opossum_spiked_file &
    #else
    #    maxJobs; nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$bamdir'/'$uid'.opossum.bam' &
    nice time python $opossum --BamFile=$bam --SoftClipsExist=True --OutFile=$opossum_normal_file
    #fi
    #echo $bamdir'/'$uid'.opossum.bam'
    else
    echo 'normal' $uid 'already done'
    fi

    platypus_spiked_file=$outDir'/'$uid'.spiked.vcf'
    platypus_normal_file=$outDir'/'$uid'.normal.vcf'

    if [ ! -f $platypus_spiked_file ]; then
        echo "Platypus running on" $bam
        #nice time python $platypus callVariants --bamFiles $bam --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $platypus_spiked_file
    else
        echo $uid 'spiked file completed'
    fi

    if [ ! -f $platypus_normal_file ]; then
        echo "Platypus running on" $bam
        #nice time python $platypus callVariants --bamFiles $matching_normal --refFile $hs37 --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o $platypus_normal_file
    else
        echo $uid 'normal file completed'
    fi

################# mutect section

split_bam="${bam%.*}".split.bam

#SplitNCigarReads
if [ ! -f $split_bam ];
then
    nice time $gatk SplitNCigarReads \
    -R $hs37 \
    -I $bam \
    -O $split_bam
fi

split_normal_bam="${matching_normal%.*}".split.bam

if [ ! -f $split_normal_bam ];
then
    nice time $gatk SplitNCigarReads \
    -R $hs37 \
    -I $bam \
    -O $split_normal_bam
fi

mutect_vcf=
mutect_normal_vcf=

$gatk Mutect2 \
-R $hs37 \
-I $split_bam \
-tumor d7bbc19b-5690-47fc-94ef-25dbade69a47 \
-pon $pon \
-O /scratch/jakutagawa/icgc/var_calls/mutect2/103c3f0e-d87e-4ae5-9727-d385c372c022.mutant.vcf
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
done

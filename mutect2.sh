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
normaldir="/scratch/jakutagawa/icgc/bams/normal"
mutantdir="/scratch/jakutagawa/icgc/bams/synthetic/spiked"
#parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/normal"
hs37="/private/groups/brookslab/reference_indices/hs37/hs37d5.fa"
hg19="/private/groups/brookslab/reference_indices/hg19/hg19.fa"
outDir="/scratch/jakutagawa/icgc/var_calls/mutect2"
pon="/scratch/jakutagawa/icgc/bams/gtex/pon/pon.vcf"
gatk="/private/groups/brookslab/jakutagawa/tools/gatk-4.0.11.0/gatk"

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
#for bam in $(find $mutantdir -mindepth 1 -name '*.with_variants.sorted.bam'); do
for bam in $(find $normaldir -mindepth 1 -name '*.v1.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    #uid=$(echo $1)
    uid=$(echo $2)
    unset IFS
    bamdir=$(dirname $bam)

    split_file=${bam%.bam}.split.bam

    if [ ! -f $split_file ]; then
        echo "SplitNCigarReads running on " $bam
        maxProcesses; nice time $gatk SplitNCigarReads -fixNDN -R $hs37 -I $bam -O $split_file &
        #echo $(find $parentdir -name $uid'.split.bam')
    else
        echo $uid' already done'
    fi

done
wait

for bam in $(find $normaldir -mindepth 1 -name '*.split.bam'); do  #'*split_picard_edited.bam'); do #'*hs37d5.bam'); do
    uid=$(basename $bam)
    IFS='.'
    set $uid
    uid=$(echo $2)
    unset IFS
    bamdir=$(dirname $bam)
    #if [ ! -d "$outDir/$uid" ]; then
        #mkdir $outDir/$uid

    echo $bam

    normal_ID_line=$(samtools view -H $bam | grep '@RG')
    IFS=$'\t'
    set $normal_ID_line
    ## !!! May need to change depending on bam header !! ##
    #GTEX bams, normal ID is at $9
#	normal_ID=$(echo $9)
    #PCAWG bams, normal ID is at $6
    normal_ID_field=$(echo $5)
    unset IFS
    normal_ID=${normal_ID_field:3}

    echo 'normal_ID ' $normal_ID
    #output_vcf=$uid.mutant.vcf
    output_vcf=$uid.normal.vcf

    #    maxProcesses; nice time gatk Mutect2 -R $hs37 -I $bamdir'/'$uid'.split.bam' -tumor $normal_ID -O $outDir'/'$uid'/'$uid'.vcf' &
    ## if running in no-match (tumor-only) mode, take out -pon argument
    if [ ! -f $outDir/$output_vcf ]; then
        if [[ $matchedOrPON = "pon"  ]]; then
            echo 'run with PON'
            echo $outDir/$output_vcf
            maxProcesses; nice time $gatk Mutect2 -R $hs37 -I $bam -tumor $normal_ID -pon $pon -O $outDir/$output_vcf &
        else
            echo 'run no PON'
            maxProcesses; nice time $gatk Mutect2 -R $hs37 -I $bam -tumor $normal_ID -O $outDir/$output_vcf &
        fi
    else
        echo $output_vcf' already exists'

    fi
    #fi
#    rm $tumordir'/'$uid'.split.bam'

    echo $uid " complete"
    done
wait

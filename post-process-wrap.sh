#!/bin/bash
#$ -cwd
####################################################################################################
# Wrapper for processing output vcfs for a given tool at a time. Finds and passes either matched 
# normal or PON file to Rscript that runs RNAedit, normal filtering and produces summary files for 
# makeCompare.R. 
# Grabs: 
#         FPKM for given DonorID 
#         Whitelist for given DonorID
#         Read coverage for tumor bam (using samtools)
#         Matched normal file names/VCF for matched mode
#    
# Input: "pon" for Panel of Normals or "" for matched normal filtering
# Output: Intermediary maf, fpkm, bed files that are used and deleted
#         Histogram plots of variants missed by tool by variant type
#         countBreakdown.txt --> contains summary numbers (filtered by editing, normal)
#            and breakdown of variants missed by tool by type
#         sensitivityFPKM.txt --> sensitivity calculated by bin range of FPKM (for makeCompare.sh)
#         sensitivityCoverage.txt --> sensitivity calculated by bin range of read coverage ^^
####################################################################################################

## PCAWG white list for all donors
#whitelist="/pod/pstore/groups/brookslab/amak/var-prediction-testing/ref-features/October_2016_whitelist_2583.snv_mnv_indel.maf"

## Subset of PCAWG whitelist for donors we have VCF output for; use "Helber code" to generate this list and subset with awk by hand
whitelist="/pod/pstore/groups/brookslab/amak/var-prediction-testing/download-lists/second-download-set-tcga/new_whitelist.maf"
#URLlist="/pod/pstore/groups/brookslab/amak/var-prediction-testing/download-lists/wgs-match-RNAseq-tumor-normal-refList-head.tsv"
tumorBamsDir="/scratch/jakutagawa/RNA-seq/realigned_bams/tumor"
URLlist="/pod/pstore/groups/brookslab/csoulette/projects/PCAWG/metadata/release_may2016.v1.4.tsv"
metadata="/pod/pstore/groups/brookslab/amak/var-prediction-testing/download-lists/second-download-set-tcga/metadata.cart.2018-05-29.json"
pcawgIDs="/pod/pstore/groups/brookslab/amak/var-prediction-testing/download-lists/second-download-set-tcga/pcawgIDs.tsv"
tophatExp="/pod/pstore/groups/brookslab/amak/var-prediction-testing/expression/star_tophat_fpkm_geneID.donor.log"
rnaEdit="/pod/pstore/groups/brookslab/amak/var-prediction-testing/ref-features/Human_AG_rnaeditsites_hg19_v2.txt"
ensemblIDs="/pod/pstore/groups/brookslab/amak/var-prediction-testing/ref-features/GRCH37_trimmed.gtf"
PON="/pod/pstore/groups/brookslab/amak/var-prediction-testing/gtex-vcfs/GTEX_PON.vcf"
## set mode to filter with matched normal or PON, options are "pon" or nothing
matchedOrPON=$1

### Uncomment out for which tool you are working with ###
## If working with VarDict ##
#parentdir="/scratch/amak/varCalls/VarDict/tumor-vcfs"
#normdir="/scratch/amak/varCalls/VarDict/normal-vcfs"
#outdir="/scratch/amak/varCalls/VarDict"
#tool="VarDict"

## If working with Platypus ##
#parentdir="/scratch/amak/varCalls/Platypus/tumor-vcfs"
#normdir="/scratch/amak/varCalls/Platypus/normal-vcfs"
#outdir="/scratch/amak/varCalls/Platypus"
#tool="Platypus"

## If working with Mutect ##
if [[ $matchedOrPON = "pon" ]]; then
    parentdir="/scratch/amak/varCalls/Mutect/PON/tumor-vcfs"
    normdir="/scratch/amak/varCalls/Mutect/normal-vcfs"
    outdir="/scratch/amak/varCalls/Mutect/PON"
    tool="Mutect2"
else
    parentdir="/scratch/amak/varCalls/Mutect/tumor-vcfs"
    normdir="/scratch/amak/varCalls/Mutect/normal-vcfs"
    outdir="/scratch/amak/varCalls/Mutect"
    tool="Mutect2"
fi
## Set out directory for plots and intermediate files
plotOuts=$(dirname $parentdir)
#For all VarDict VCFs --> will find all corresponding VCFs from other tools
for vcf in $(find $parentdir -name '*.vcf'); do
    echo 
    uid=$(basename $vcf)
    IFS='_.'; set $uid; uid=$(echo $1) 
    IFS=''

    echo 'uid' $uid
    tumorBam=$(find $tumorBamsDir -name "*$uid*.bam")
    echo 'tumorBam'  $tumorBam
    # Grab the case ID from the metadata file and use it to find the match normal file name. Use the case IDs to find the line from the URL list containing the Donor ID and match-normal ID
    caseID=$(grep -A 10 $uid $metadata | grep case_id)

    IFS=':'; set $caseID; caseID=$(echo $2 | sed 's/", //g' | sed 's/ "//g')
    IFS=''
    echo 'Case ID: ' $caseID
    # Grab filename for match normal
    normalID=$(grep -B 10 $caseID $metadata | grep submitter_id | grep -v $uid)  
    IFS=':'; set $normalID; normalID=$(echo $2 | sed 's/", //g' | sed 's/ "//g')
    echo "NORMAL ID:" $normalID

    donorLine=$(grep $caseID $URLlist)
#    echo $donorLine



    #Grab Donor ID

    donorID=$(echo $donorLine | awk 'BEGIN {FS="\t"}; {print $5}')
    echo 'Donor ID: ' $donorID

    ## Helber code for making a donor-list for filtering, makes list of donors so can filter to reduce dimensions of PCAWG whitelist ##

#    if grep -q $donorID '/pod/pstore/groups/brookslab/amak/scripts/Variant-Calling/donor-list.txt'; then
#	continue
#    else
#	echo $donorID >> 'donor-list.txt'
#    fi

    #Write whitelist MAF to temporary file 
    donorMaf=$parentdir'/'$donorID'.maf'

    time grep $donorID $whitelist > $parentdir'/'$donorID'.maf'


    #Grab expression data for Donor
    time awk -v col=$donorID 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $1, "\t", $2, "\t", $c}' $tophatExp > $parentdir'/'$donorID'.fpkm'
    sed -i 1d $parentdir'/'$donorID'.fpkm'
    ## code for when file id is the same as UID (ie for PCAWG bams); current code is for reheaded bams 
    #Grab Match Normal
#    STARNormal=$(echo $donorLine | awk 'BEGIN {FS="\t"}; {print $73}')
#    STARTumor=$(echo $donorLine | awk 'BEGIN {FS="\t"}; {print $87}')
#    tophatNormal=$(echo $donorLine | awk 'BEGIN {FS="\t"}; {print $77}')
#    tophatTumor=$(echo $donorLine | awk 'BEGIN {FS="\t"}; {print $91}')
#    echo $tophatNormal
#    IFS=.
#    set $STARNormal
#    normal=$(echo $2 | sed 's/", //g' | sed 's/ "//g')   
#    set $STARTumor
#    tumor=$(echo $2 | sed 's/", //g' | sed 's/ "//g')
#    echo 'tumor' $tumor
#    unset IFS
#    findStar='*'$STARNormal'*'
    if [[ $matchedOrPON = "pon" ]]; then 
	echo "in PON mode!"
	VCFN=$PON

	ponMode=true
    else
	echo "in matched normal mode!"
	find $normdir -mindepth 1 -type f -name "*${normalID}*.vcf" -print0 | xargs -0 grep -v "\#" > $parentdir'/'$donorID'.normal.tmp'
	VCFN="$parentdir/$donorID.normal.tmp"
	ponMode=false
    fi
    grep -v "\#" $vcf > $parentdir'/'$donorID'.tumor.tmp'
    VCFT="$parentdir/$donorID.tumor.tmp"


    echo 'donor MAF: ' $donorMaf
    echo 'donor id ' $donorID

    ## R script takes temp maf and temp FPKM for given donor and writes 
    # countBreakdown.txt <- counts of variants missed compared to whitelist
    # sensitivityFPKM.txt <- sensitivity vs FPKM data frame
    # found/notFound.bed <- bed files for extracting read coverage at variants found/not found
    time Rscript compareWhitelist.R $parentdir $parentdir'/'$donorID'.maf' $VCFT $VCFN $rnaEdit $parentdir'/'$donorID'.fpkm' $ensemblIDs $donorID $outdir $ponMode >> post-process.log
#    echo "time Rscript compareWhitelist.R $parentdir $parentdir/$donorID.maf $VCFT $VCFN $rnaEdit $parentdir/$donorID.fpkm $donorID"
    
    echo 
    ## call bedtools coverage to get coverage 
    if [[ ! -f $outdir'/found.cov.txt' ]]; then
	samtools bedcov $parentdir'/'$donorID'.found.bed' $tumorBam > $parentdir'/'$donorID'.found.cov.txt'
	samtools bedcov $parentdir'/'$donorID'.notFound.bed' $tumorBam > $parentdir'/'$donorID'.notFound.cov.txt'
	#    bedtools coverage -abam $tumorBam -b $parentdir'/'$donorID'.found.bed' > $parentdir'/'$donorID'.found.cov.txt' -counts
    fi
#    bedtools coverage -abam $tumorBam -b $parentdir'/'$donorID'.notFound.bed' > $parentdir'/'$donorID'.notFound.cov.txt' -counts
    ## Wait until R script is done and delete temporary donor maf and fpkm files ##
    wait

    find $parentdir -name "*.tmp" -type f -delete
    find $parentdir -name "*.maf" -type f -delete
    find $parentdir -name "*.fpkm" -type f -delete
#    find $parentdir -name "*.bed" -type f -delete

 


done




## Add header to file outlining breakdown of missed variants compared to whitelist ## 
echo -e "Donor ID \t Sensitivity \t numFiltGermline \t numFiltEdits \t 3'UTR \t 5'UTR \t RNA \t lincRNA \t Silent \t Splice_Site \t De_novo_Start_InFrame \t De_novo_Start_OutOfFrame \t Missense_Mutation \t Nonsense_Mutation \t Nonstop_Mutation \t In_Frame_Del \t In_Frame_Ins \t Frame_Shift_Del \t Frame_Shift_Ins \t Start_Codon_Del \t Start_Codon_Ins \t Start_Codon_SNP \t Stop_Codon_Del \t Stop_Codon_Ins \t 5'Flank \t Intron \t IGR" | cat - $outdir'/countBreakdown.txt' > tmp && mv tmp $outdir'/countBreakdown.txt' 
echo -e "Donor ID \t Total \t 3'UTR \t 5'UTR \t RNA \t lincRNA \t Silent \t Splice_Site \t De_novo_Start_InFrame \t De_novo_Start_OutOfFrame \t Missense_Mutation \t Nonsense_Mutation \t Nonstop_Mutation \t In_Frame_Del \t In_Frame_Ins \t Frame_Shift_Del \t Frame_Shift_Ins \t Start_Codon_Del \t Start_Codon_Ins \t Start_Codon_SNP \t Stop_Codon_Del \t Stop_Codon_Ins \t 5'Flank \t Intron \t IGR" | cat - $outdir'/totalBreakdown.txt' > tmp && mv tmp $outdir'/totalBreakdown.txt'

## Plot count breakdown and sensitivity vs FPKM 
Rscript histogramBreakdown.R $outdir $tool

## put all found coverage files together to pass to makeCompiled.R and delete intermediary files
if [[ ! -f $outdir'/found.cov.txt' ]]; then
    find $parentdir -type f -name '*.found.cov.txt' | xargs cat >> $outdir'/found.cov.txt'
    find $parentdir -type f -name '*.notFound.cov.txt' | xargs cat >> $outdir'/notFound.cov.txt'
    find $parentdir -type f -name '*.found.cov.txt' -delete
    find $parentdir -type f -name '*.notFound.cov.txt' -delete
fi

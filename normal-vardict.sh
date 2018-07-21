#!/bin/bash
#$ -cwd 


#outputs concatenated vardict file for now, hopefully vcf soon

#parentdir=$1 #"/scratch/jakutagawa/RNA-seq/tumor_bams/reindex/"
parentdir="/scratch/jakutagawa/RNA-seq/realigned_bams/normal"
hg19="/pod/pstore/groups/brookslab/reference_indices/hg19/hg19.fa"
hs37="/pod/pstore/groups/brookslab/reference_indices/hs37/hs37d5.fa"
cdnahg19="/pod/pstore/groups/brookslab/reference_indices/hg19/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa"
splitbeds="/scratch/amak/varCalls/VarDict/wg-beds"
outDir="/scratch/amak/varCalls/VarDict/normal-vcfs"
bedList=()

function maxJobs {
    # Waits until there are less than 'numJobs' jobs running before starting a new job        
    while [ $(jobs | wc -l) -ge 8 ]; do
        echo 'waiting'
        sleep 5
    done
}


for bed in $(ls $splitbeds); do 

    bedList+=($bed)
done
echo ${bedList[@]}

for bam in $(find $parentdir -mindepth 1 -name '*hs37d.bam'); do
    uid=$(basename $bam)

    mkdir $outDir/$uid
#    cd $outDir/$uid
    echo "VarDict running on " $uid
    for bed in ${bedList[@]}; do
#	echo $bed
	
	maxJobs; nice time vardict -D -G $hs37 -f 0.01 -N $uid -b $bam -c 1 -S 2 -E 3 $splitbeds/$bed >> $outDir/$uid/$uid".out" &

    done
    
    python vardict2vcf.py $outDir/$uid/$uid".out" > $outDir/$uid/$uid".vcf"
    echo $uid " complete"
#    cat $outDir/$uid/* > $outDir/$uid/
done
wait

#time /scratch/amak/packages/VarDictJava/VarDict/vardict -D -G /pod/pstore/groups/brookslab/reference_indices/hg19/hg19.fa -f 0.01 -N PCAWG.ba92c434-3604-4b85-bc76-3bbe5c44253f.TopHat2.v1.chr.sorted -b /scratch/jakutagawa/RNA-seq/tumor_bams/reindex/c1b68b54-0258-470b-b6e2-b3f558bb1293/PCAWG.ba92c434-3604-4b85-bc76-3bbe5c44253f.TopHat2.v1.chr.sorted -c 1 -S 2 -E 3 /scratch/jakutagawa/RNA-seq/tumor_bams/reindex/hg19_5k_150bpOL_seg.txt.2 | /scratch/amak/packages/VarDictJava/VarDict/teststrandbias.R | /scratch/amak/packages/VarDictJava/VarDict/var2vcf_valid.pl -N PCAWG.ba92c434-3604-4b85-bc76-3bbe5c44253f.TopHat2.v1.chr.sorted -E -f 0.01 > /scratch/amak/varCalls/VarDictJava/c1b68b54-0258-470b-b6e2-b3f558bb1293/PCAWG.ba92c434-3604-4b85-bc76-3bbe5c44253f.TopHat2.v1.chr.sorted.2.vcf &



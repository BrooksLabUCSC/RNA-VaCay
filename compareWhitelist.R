#!/usr/bin/env Rscript

library("dplyr")
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Supply working directory, whitelist, tumor, normal VCFs, and .n", call.=FALSE)
}


setwd(args[1])
parentdir=args[1]
whitelist=read.delim(args[2], sep="\t", header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
tumorVCF=read.table(args[3], sep='\t', header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
normalVCF=read.table(args[4], sep='\t', header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
rnaEdit=read.table(args[5], sep='\t', header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE)
doFPKM=read.table(args[6], sep='\t', header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
ensembl=read.table(args[7], sep='\t', header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
donorID=args[8]
outdir=args[9]
ponMode=args[10]

print(donorID)
print("tumor")
print(args[3])
print("normal")
print(args[4])

########################## Check which tool VCF came from ############################

if (grepl("Platypus", args[1])==TRUE) {
   tool <- "Platypus"
   ############# Add 'chr' to column 1 ###############
   normalVCF$V1 <- paste('chr', normalVCF$V1, sep="")
   tumorVCF$V1 <- paste('chr', tumorVCF$V1, sep="")   
   print("In PLATYPUS")
   library("tidyr")
   print(head(tumorVCF))
   print(head(normalVCF))
   # Separate 'INFO' column and extract VAF and Coverage #

#   normalVCF <- normalVCF %>% extract(V10, c("one", "VC"), regex="(.*):(.*)$")
#   tumorVCF <- tumorVCF %>% extract(V10, c("one", "VC"), regex="(.*):(.*)$")

#   normalVCF <- normalVCF %>% separate(V8, c("V8", "VAF", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "TC", "TCF", "TCR"), ";")
#   tumorVCF <- tumorVCF %>% separate(V8, c("V8", "VAF", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "TC", "TCF", "TCR"), ";")

#   normalVCF <- select(normalVCF, -8, -(10:21), -25, -one)
#   tumorVCF <- select(tumorVCF, -8, -(10:21), -25, -one)

    normalVCF <- select(normalVCF, (1:8))
   tumorVCF <- select(tumorVCF, (1:8))


   ########## Remove unnecessary characters ###########
#   tumorVCF$VAF <- gsub("[TCFR\\=]", "", tumorVCF$VAF)
#   tumorVCF$TC <- gsub("[TCFR\\=]", "", tumorVCF$TC)
#   tumorVCF$TCF <- gsub("[TCFR\\=]", "", tumorVCF$TCF)
#   tumorVCF$TCR <- gsub("[TCFR\\=]", "", tumorVCF$TCR)

#   normalVCF$VAF <- gsub("[TCFR\\=]", "", normalVCF$VAF)
#   normalVCF$TC <- gsub("[TCFR\\=]", "", normalVCF$TC)
#   normalVCF$TCF <- gsub("[TCFR\\=]", "", normalVCF$TCF)
#   normalVCF$TCR <- gsub("[TCFR\\=]", "", normalVCF$TCR)
   print("tumor VCF Platypus")
   print(head(tumorVCF))
   normalVCF <- as.data.frame(normalVCF)
   tumorVCF <- as.data.frame(tumorVCF)
   ############### Change column names ################
   colnames(tumorVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF") #, "Total.Coverage", "Forward.Coverage", "Reverse.Coverage", "Alt.Coverage")
   colnames(normalVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF") #, "Total.Coverage", "Forward.Coverage", "Reverse.Coverage", "Alt.Coverage")

} else if (grepl("VarDict", args[1])==TRUE) {
  tool <- "VarDict"
  ############### Change column names ################
  normalVCF$V1 <- paste('chr', normalVCF$V1, sep="")	
  tumorVCF$V1 <- paste('chr', tumorVCF$V1, sep="")

  library("tidyr")
  normalVCF <- normalVCF %>% separate(V8, c("VAF", "TC", "VC"), ";")
  tumorVCF <- tumorVCF %>% separate(V8, c("VAF", "TC", "VC"), ";")
  tumorVCF$VAF <- gsub("[VATCFR\\=]", "", tumorVCF$VAF)
  tumorVCF$TC <- gsub("[VATCFR\\=]", "", tumorVCF$TC)
  tumorVCF$VC <- gsub("[VATCFR\\=]", "", tumorVCF$VC)

  normalVCF$VAF <- gsub("[VATCFR\\=]", "", normalVCF$VAF)
  normalVCF$TC <- gsub("[VATCFR\\=]", "", normalVCF$TC)
  normalVCF$VC <- gsub("[VATCFR\\=]", "", normalVCF$VC)

  normalVCF <- as.data.frame(normalVCF)
  tumorVCF <- as.data.frame(tumorVCF)
  print("tumor VCF VarDict")
  print(head(tumorVCF))
  print(head(normalVCF))
  normalVCF <- select(normalVCF, -11, -12)
  colnames(tumorVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF", "Total.Coverage", "Alt.Coverage")  
  colnames(normalVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF", "Total.Coverage", "Alt.Coverage")
  print(head(normalVCF))

  ####### Change varDict to match whitelist variant type ###########
#  whitelist <- replace(whitelist$V7, whitelist$V7=="SNP","SNV")
#  whitelist <- replace(whitelist$V7, whitelist$V7=="INS","Insertion")
#  whitelist <- replace(whitelist$V7, whitelist$V7=="DEL","Deletion")
#  whitelist <- replace(whitelist$V7, whitelist$V7=="DNP","Complex")
} else {
  tool <- "Mutect2"
  print("in Mutect")
  ############### Change column names ################
  normalVCF$V1 <- paste('chr', normalVCF$V1, sep="")
  tumorVCF$V1 <- paste('chr', tumorVCF$V1, sep="")

  library("tidyr")
  normalVCF <- normalVCF %>% separate(V8, c("VAF", "TC", "VC"), ";")
  tumorVCF <- tumorVCF %>% separate(V8, c("VAF", "TC", "VC"), ";")
  tumorVCF$VAF <- gsub("[VATCFR\\=]", "", tumorVCF$VAF)
  tumorVCF$TC <- gsub("[VATCFR\\=]", "", tumorVCF$TC)
  tumorVCF$VC <- gsub("[VATCFR\\=]", "", tumorVCF$VC)

  normalVCF$VAF <- gsub("[VATCFR\\=]", "", normalVCF$VAF)
  normalVCF$TC <- gsub("[VATCFR\\=]", "", normalVCF$TC)
  normalVCF$VC <- gsub("[VATCFR\\=]", "", normalVCF$VC)

  normalVCF <- as.data.frame(normalVCF)
  tumorVCF <- as.data.frame(tumorVCF)
  print("tumor VCF Mutect")
  print(head(tumorVCF))
  print(head(normalVCF))
  normalVCF <- select(normalVCF, -11, -12)
  tumorVCF  <- select(tumorVCF, -11, -12)
  colnames(tumorVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF", "Total.Coverage", "Alt.Coverage")
  colnames(normalVCF) <- c("Chr", "Start", "ID", "Ref", "Alt", "Quality", "Filter", "VAF", "Total.Coverage", "Alt.Coverage")
  print('edited colnames')
  print(head(normalVCF))
  print(head(tumorVCF))
}


########################## Filter out normals ########################## 
tumorVCF$Start <- as.numeric(tumorVCF$Start)
normalVCF$Start <- as.numeric(normalVCF$Start)
#Filter out incomplete lines

normalVCF <- normalVCF %>% filter(!is.na(normalVCF$Start))
tumorVCF <- tumorVCF %>% filter(!is.na(tumorVCF$Start))

if (isTRUE(all.equal(ponMode, 'true')) & isTRUE(all.equal(tool, 'Mutect2'))) {
   print("Mutect2 no PON")
   numFiltGermline=0
   somatic <- tumorVCF
} else {
  print("Filtering with PON")
  # Filter out normal calls
  somatic <- tumorVCF %>% anti_join(normalVCF, by="Start", "Chr")
  numFiltGermline <- dim(tumorVCF)[1] - dim(somatic)[1]
 
}
########################## Filter out RNA editing sites ##########################
rnaEdit$position <- as.numeric(rnaEdit$position)

somaticMinusEdits <- somatic %>% anti_join(rnaEdit, by=c("Chr" = "chromosome", "Start" = "position") )
head(somaticMinusEdits)
numFiltEdits <- dim(somatic)[1] - dim(somaticMinusEdits)[1]
numFiltEdits

########################## Filter out low mapQ scores ##########################
somaticMQ50 <- filter(somaticMinusEdits, somaticMinusEdits$Quality >= 30)
numFiltMQ50 <- dim(somaticMinusEdits)[1] - dim(somaticMQ50)[1]

########################## Compare to whitelist MAF ##########################
# V3 is start position, V6 variant location, V7 variant type (SNP), V17 is probably VAF?, V42 cancer


whitelist$V2 <- paste('chr', whitelist$V2, sep="")

foundInMaf <- somaticMinusEdits %>% semi_join(whitelist, by=c("Chr"="V2", "Start"="V3"))
print("somaticMinusEdits:")
dim(somaticMinusEdits)
dim(foundInMaf)
#head(foundInMaf)
#correctClass <- somaticMinusEdits %>% semi_join(whitelist, by=c("Chr"="V2", "Start"="V3", "Variant"="V7"))

totalMaf <- dim(whitelist)[1]
numInMaf <- dim(foundInMaf)[1]
#sensitivity <- numInMaf / (numInMaf + totalMaf)
#####specificity <- 

#numCorrectClass <- dim(correctClass)[1]
#sensitivityCorrectClass <- numCorrectClass / (numCorrectClass + totalMaf)

### Whitelist annotations for variants found by caller ###
Maflocs <- whitelist %>% semi_join(somaticMinusEdits, by=c("V2"="Chr", "V3"="Start"))
foundBreakdown <- Maflocs %>% group_by(V6) %>% summarize(count=n())
rownames(foundBreakdown) <- foundBreakdown$V6
print('foundBreakdown')
print(as.data.frame(foundBreakdown))

### whitelist annotations for variants not found by caller ###
MafNotFound <- whitelist %>% anti_join(somaticMinusEdits, by=c("V2"="Chr", "V3"="Start"))
MafNotFound <- as.data.frame(MafNotFound)

breakdown <- as.data.frame(MafNotFound %>% group_by(V6) %>% summarise(count=n()))
rownames(breakdown) <- breakdown$V6
print('breakdown')
print(as.data.frame(breakdown))

# Filter out intronic, intergenic
notFoundMinusIntron <- MafNotFound %>% filter(V6!="Intron", V6!="IGR")

# Get sensitivity minus intronic/intergenic
MafMinusIGR <- dim(notFoundMinusIntron)[1]
sensitivity <- numInMaf / (numInMaf + MafMinusIGR)

## Print sanity check to R.log ##
print("numInMaf")
print(numInMaf)
print("totalMAF")
print(MafMinusIGR)
print("sensitivity")
print(sensitivity)

############## Write sensitivity, specificity, breakdown, filtered out to file ##############

outList = list(donorID, sensitivity, numFiltGermline, numFiltEdits, 
	as.numeric(breakdown["3'UTR",][2]), 	      
	as.numeric(breakdown["5'UTR",][2]), 
	as.numeric(breakdown["RNA",][2]), 
	as.numeric(breakdown["lincRNA",][2]), 
	as.numeric(breakdown["Silent",][2]), 
	as.numeric(breakdown["Splice_Site",][2]),
	as.numeric(breakdown["De_novo_Start_InFrame",][2]),
	as.numeric(breakdown["De_novo_Start_OutOfFrame",][2]),
	as.numeric(breakdown["Missense_Mutation",][2]), 
	as.numeric(breakdown["Nonsense_Mutation",][2]), 
	as.numeric(breakdown["Nonstop_Mutation",][2]), 
	as.numeric(breakdown["In_Frame_Del",][2]),
	as.numeric(breakdown["In_Frame_Ins",][2]),
	as.numeric(breakdown["Frame_Shift_Del",][2]),
	as.numeric(breakdown["Frame_Shift_Ins",][2]),
	as.numeric(breakdown["Start_Codon_Del",][2]),
	as.numeric(breakdown["Start_Codon_Ins",][2]),
	as.numeric(breakdown["Start_Codon_SNP",][2]),
	as.numeric(breakdown["Stop_Codon_Del",][2]),
	as.numeric(breakdown["Stop_Codon_Ins",][2]),
	as.numeric(breakdown["5'Flank",][2]),
	as.numeric(breakdown["Intron",][2]), 
	as.numeric(breakdown["IGR",][2]))

countsOut=paste(outdir, "/countBreakdown.txt", sep="")
write.table(outList, file=countsOut, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)

 
########### Get sensitivity per Variant Type ###################
variant.types <- list("Donor ID", "Sensitivity", "3'UTR", "5'UTR", "RNA", "lincRNA", "Silent", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins", "5'Flank", "Intron", "IGR")


sensList <- list()
i <- 1
for (var in variant.types) {
    print(var)
    sensitivity <- as.numeric(foundBreakdown[var,][2]) / (as.numeric(foundBreakdown[var,][2]) + as.numeric(breakdown[var,][2]))
    sensList[[i]] <- sensitivity
    print(sensitivity)
    i <- i + 1

}
sensList[1] <- donorID
sensList[2] <- sensitivity

sensOut=paste(outdir, "/sensitivityBreakdown.txt", sep="")
write.table(sensList, file=sensOut, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)

## make df with totals of each variant type for a given donor
totList <- list()
i <- 1
for (var in variant.types) {
    print(var)
    print('total')
    print(breakdown[var,][2])
    found <- as.numeric(foundBreakdown[var,][2])
    notFound <- as.numeric(breakdown[var,][2])
    found <- replace(found, is.na(found), 0)
    notFound <- replace(notFound, is.na(notFound), 0)
    total <- found + notFound
    totList[[i]] <- total
    print(total)
    i <- i + 1

}
totList[1] <- donorID
totList[2] <- total

totOut=paste(outdir, "/totalBreakdown.txt", sep="")
write.table(totList, file=totOut, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)


########################## Compare FPKM ##########################
## NOTE: tophat/STAR FPKM only has gene-level positions (must be average FPKM across gene) so introns are picked up ##

#### set colnames for donor FPKM and separate out EnsemblID ####
colnames(doFPKM) <- c("Ensembl", "GeneID", "FPKM")
doFPKM <- doFPKM %>% separate(Ensembl, c("EnsemblID", "Version"), sep="\\.")
doFPKM <- doFPKM %>% select(-Version)
#### set colnames for ensemblID position df and add 'chr' to chrom ####
colnames(ensembl) <- c("Chr", "Type", "Start", "End", "EnsemblID")
ensembl$Chr <- paste('chr', ensembl$Chr, sep="")

#### join start/end positions to donor FPKM by ensembl ID ####
ensFPKM <- doFPKM %>% inner_join(ensembl, by="EnsemblID")
#doFPKM <- doFPKM[-1,]

#### clean up whitelist ####
doMAF <- whitelist %>% select(1:10, 17, 29, 42, 43)

#### join donor MAF by start in range of start/end in donor MAF ####

## remove variants in introns and intergenic regions (increases speed) ##
doMAF <- doMAF %>% filter(V6 != "5'Flank" & V6 != "Intron" & V6 != "IGR")

# this way is slower but equivalent #
#mafFPKM <- doMAF %>% mutate(dummy="TRUE") %>% left_join(ensFPKM %>% mutate(dummy="TRUE")) %>% filter(V2==Chr, V3>=Start, V3<=End) %>% select(-dummy)


## This removes 3891 variant sites from whitelist, likely not expressed
mafFPKM <- doMAF %>% rowwise() %>% mutate(s=list(ensFPKM %>% filter(V2==Chr, V3>Start, V3<End))) %>% tidyr::unnest()
mafFPKM <- as.data.frame(mafFPKM %>% distinct(mafFPKM, V3, .keep_all=TRUE))


#### log transform Maf/FPKM combined dataframe ####
# if FPKM is 0
mafFPKM$logFPKM[mafFPKM$FPKM == 0] = log2(1)
# if FPKM is positive, take log2
mafFPKM$logFPKM[mafFPKM$FPKM > 0] = log2(mafFPKM$FPKM[mafFPKM$FPKM >0])
# if FPKM is negative, take abs(log2) and put negative back on
mafFPKM$logFPKM[mafFPKM$FPKM < 0] = -log2(abs(mafFPKM$FPKM[mafFPKM$FPKM <0]))


#### Set bins by FPKM value on combined MAF/FPKM df ####
#mafFPKM <- mafFPKM %>% mutate(bin=ntile(as.numeric(FPKM), 10))
mafFPKM$bin <- cut(mafFPKM$logFPKM, breaks=c(-12.5, -12, -11.5, -11, -10.5, -10, -9.5, -9, -8.5, -8, -7.5, -7, -6.5, -6, -5.5, -5, -4.5, -4,-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5), labels=c("-12.5", "-12", "-11.5", "-11", "-10.5", "-10", "-9.5", "-9", "-8.5", "-8", "-7.5", "-7", "-6.5", "-6", "-5.5", "-5", "-4.5", "-4", "-3.5", "-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", "12", "12.5"))

#### Reselect foundInMAF and get breakdowns for found and total ####

foundMaf <- doMAF %>% semi_join(somaticMinusEdits, by=c("V2"="Chr", "V3"="Start"))
notfoundMaf <- doMAF %>% anti_join(somaticMinusEdits, by=c("V2"="Chr", "V3"="Start"))

foundMafFPKM <- mafFPKM %>% semi_join(somaticMinusEdits, by=c("V2"="Chr", "V3"="Start"))
## Makes column for counts
foundMafFPKM <- foundMafFPKM %>% group_by(bin) %>% mutate(count=n())
print("foundMaf")
print(head(foundMaf))

## Get counts by FPKM bin for not found (sensFPKM) and all (foundsens) variants ##
sensFPKM <- as.data.frame(mafFPKM %>% group_by(bin) %>% summarise(total=n()))
foundsens <- as.data.frame(foundMafFPKM %>% group_by(bin) %>% summarise(found=n()))

## Check means for each bin to see if associating with correct bin ##
meanFPKM <- as.data.frame(mafFPKM %>% group_by(bin) %>% summarise(mean=mean(logFPKM)))

## Join counts of found and not found variants by FPKM bin ##
sensFPKM <- sensFPKM %>% left_join(foundsens, by="bin")

print("sensFPKM")
print(head(sensFPKM))
## Add sensitivity column and replace NAs with 0s ##
sensFPKM <- as.data.frame(sensFPKM %>% group_by(bin) %>% mutate(sensitivity=found / (total)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))

# paste outdir to appropriate file name and append FPKM with counts in each bin to file
# this will be used later in makeCompiled.R to make sensitivity vs FPKM for all samples
sensitivityFPKM=paste(outdir, "/sensitivityFPKM.txt", sep="")
write.table(sensFPKM, file=sensitivityFPKM, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)

##################### Append to bed files for found and not found to extract coverage in makeCompiled #################

#make bed of found
foundBed <- as.data.frame(foundMaf %>% select(V2, V3, V4, V1, V29, V5))
notFoundBed <- as.data.frame(notfoundMaf %>% select(V2, V3, V4, V1, V29, V5))
foundBed$V2 <- gsub("[chr\\=]", "", foundBed$V2)
notFoundBed$V2 <- gsub("[chr\\=]", "", notFoundBed$V2)
##### add 1 to all end positions if start and end are equal #####
# making end position start+1 because otherwise samtools gives sum of coverage at each nuc
foundBed[,3] <- foundBed[,2] + 1
notFoundBed[,3] <- notFoundBed[,2] + 1


foundBedFile=paste(parentdir, "/", donorID, ".found.bed", sep="")
notFoundBedFile=paste(parentdir, "/", donorID, ".notFound.bed", sep="")
print(head(foundBed))
print("foundBedFile")
print(foundBedFile)
write.table(foundBed, file=foundBedFile, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)
write.table(notFoundBed, file=notFoundBedFile, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)




#### This is temporary, removed 4229 genes from MAF ####
mafFPKM <- doMAF %>% semi_join(doFPKM, by=c("V1"="GeneID"))
mafFPKM <- mafFPKM %>% semi_join(doFPKM, by=c("V1"="GeneID"))
#filt.DO38 <- DO38901 %>% filter(DO38901, DO38901 > 1)
#filtMafExp <- filt.DO38 %>% semi_join(whiteMaf, by=c("genes"="Hugo_Symbol"))
#filtMaf <- whiteMaf %>% semi_join(filt.DO38, by=c("Hugo_Symbol"="genes"))

#filtMaf %>% group_by(Variant_Classification) %>% summarize(count=n())

#newPlat %>% mutate(bin=ntile(as.numeric(V22), 3))

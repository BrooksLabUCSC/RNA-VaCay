#!/usr/bin/env Rscript

library("dplyr")
library(ggplot2)
library(reshape2)
setwd("/scratch/amak/varCalls")

# Change these directories to fit your file system to where outDir was for post-process-wrap.sh
platDir <- "/scratch/amak/varCalls/Platypus/"
varDir <- "/scratch/amak/varCalls/VarDict/"
mutDir <- "/scratch/amak/varCalls/Mutect/"
#

platCounts <- read.delim(paste(platDir, 'countBreakdown.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
varCounts <- read.delim(paste(varDir, 'countBreakdown.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
mutCounts <- read.delim(paste(mutDir, 'countBreakdown.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

# total number of variants in whitelist
varTotal <- read.delim(paste(varDir, 'totalBreakdown.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

platSens <- read.delim(paste(platDir, 'sensitivityFPKM.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
varSens <- read.delim(paste(varDir, 'sensitivityFPKM.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
mutSens <- read.delim(paste(mutDir, 'sensitivityFPKM.txt', sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

platCovFound <- read.delim(paste(platDir, 'found.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
platCovNotFound <- read.delim(paste(platDir, 'notFound.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)

varCovFound <- read.delim(paste(varDir, 'found.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
varCovNotFound <- read.delim(paste(varDir, 'notFound.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)

mutCovFound <- read.delim(paste(mutDir, 'found.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
mutCovNotFound <- read.delim(paste(mutDir, 'notFound.cov.txt', sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)

## Process for Platypus ##
platCounts[is.na(platCounts)] <- 0
platSensitivity <- as.numeric(platCounts %>% summarize(Mean=mean(Sensitivity)))
platFiltGermline <- as.numeric(platCounts %>% summarize(Mean=mean(numFiltGermline)))
platFiltEdits <- as.numeric(platCounts %>% summarize(Mean=mean(numFiltEdits)))
# add label column for tool type so can concat later
platCovFound <- platCovFound %>% mutate(Tool="Platypus")
platCovNotFound <- platCovNotFound %>% mutate(Tool="Platypus")
## Print overall stats ##
print("Platypus Average Sensitivity")
print(platSensitivity)
print("Platypus Average Filtered PON")
print(platFiltGermline)
print("Platypus Filtered RNAEdits")
print(platFiltEdits)

## Process for VarDict ##
varCounts[is.na(varCounts)] <- 0
varSensitivity <- as.numeric(varCounts %>% summarize(Mean=mean(Sensitivity)))
varFiltGermline <- as.numeric(varCounts %>% summarize(Mean=mean(numFiltGermline)))
varFiltEdits <- as.numeric(varCounts %>% summarize(Mean=mean(numFiltEdits)))
varCovFound <-	varCovFound %>% mutate(Tool="VarDict")
varCovNotFound	<- varCovNotFound %>% mutate(Tool="VarDict")
print("VarDict Average Sensitivity")
print(varSensitivity)
print("VarDict Average Filtered PON")
print(varFiltGermline)
print("VarDict Average Filtered RNAEdits")
print(varFiltEdits)

## Process for Mutect ##
mutCounts[is.na(mutCounts)] <- 0
mutSensitivity <- as.numeric(mutCounts %>% summarize(Mean=mean(Sensitivity)))
mutFiltGermline <- as.numeric(mutCounts %>% summarize(Mean=mean(numFiltGermline)))
mutFiltEdits <- as.numeric(mutCounts %>% summarize(Mean=mean(numFiltEdits)))
mutCovFound <- mutCovFound %>% mutate(Tool="Mutect2")
mutCovNotFound	 <- mutCovNotFound %>% mutate(Tool="Mutect2")
print("Mutect Average Sensitivity")
print(mutSensitivity)
print("Mutect Average Filtered PON")
print(mutFiltGermline)
print("Mutect Average Filtered RNAEdits")
print(mutFiltEdits)


## Remove columns and add label column to each dataset ##
platCounts <- select(platCounts, -Donor.ID, -Sensitivity, -numFiltGermline, -numFiltEdits) #%>% mutate(Tool="Platypus")
varCounts <- select(varCounts, -Donor.ID, -Sensitivity, -numFiltGermline, -numFiltEdits) #%>% mutate(Tool="VarDict")
mutCounts <- select(mutCounts, -Donor.ID, -Sensitivity, -numFiltGermline, -numFiltEdits) #%>% mutate(Tool="Mutect2")

## sum down each column(variant type) to get total not found for all bams
platCounts <- platCounts %>% summarize_all(sum) 
varCounts <- varCounts %>% summarize_all(sum) 
mutCounts <- mutCounts %>% summarize_all(sum)
varTotal <- varTotal %>% mutate_all(funs(ifelse(is.na(.), 0, .))) %>% select(-Donor.ID, -Total)
print("varCounts")
print(head(varCounts))
## sum down each column (variant type) to get total variants in whitelist for all bams
varTotal <- varTotal %>% summarize_all(sum)
print("varTotal")
print(varTotal)
library(tidyr)
## sort by lowest to highest sum
platCounts <- platCounts %>% gather(Variant.Type, Sum) %>% arrange(desc(Sum))
varCounts <- varCounts %>% gather(Variant.Type, Sum) %>% arrange(desc(Sum))
mutCounts <- mutCounts %>% gather(Variant.Type, Sum) %>% arrange(desc(Sum))
varTotal <- varTotal %>% gather(Variant.Type, Total) %>% arrange(desc(Total))
positions <- varCounts$Variant.Type

## to print total num of variants in same order as varTotal and add column to retain tool type before bind
platCounts <- platCounts %>% full_join(varTotal, by="Variant.Type") %>% mutate(Tool="Platypus")
varCounts <- varCounts %>% full_join(varTotal, by="Variant.Type") %>% mutate(Tool="VarDict")
mutCounts <- mutCounts %>% full_join(varTotal, by="Variant.Type") %>% mutate(Tool="Mutect2")

## Bind data frames together (append) ##
sumBreakdown <- bind_rows(platCounts, varCounts, mutCounts)
#sumBreakdown$Tool <- factor(sumBreakdown$Tool, levels = c("VarDict", "Platypus", "Mutect2"))
maxSum=max(sumBreakdown$Sum)
print(sumBreakdown)
######## Plot sum with IGR, 5'Flank, and Introns ##########
withIGRHisto <- ggplot(sumBreakdown, aes(x=Variant.Type, y=Sum, fill=Tool)) + geom_bar(stat="identity", position="dodge") + coord_flip() + geom_text(aes(label=paste(Sum, " of ", Total)), hjust=-0.3, size=3, position = position_dodge(width = 1)) + theme_classic() + guides(fill = guide_legend(reverse = TRUE)) + scale_x_discrete(limits=positions) + scale_y_continuous("Count", expand = c(0, 0), limits=c(0,maxSum+6000))
ggsave(withIGRHisto, filename="HistoWithIGR_4.jpg")

########## Plot Sum without IGR, 5'Flank, and Introns ###########
sumBreakdown <- sumBreakdown %>% filter(Variant.Type != "IGR" & Variant.Type != "Intron" & Variant.Type !="X5.Flank")
print(sumBreakdown)
#sumBreakdown <- sumBreakdown %>% gather(Variant.Type, Mean) %>% arrange(desc(Sum))
print(positions)
positions<- positions[!is.element(positions, c("IGR", "Intron", "X5.Flank"))] #%>% filter(Variant.Type != "IGR" & Variant.Type != "Intron" & Variant.Type != "X5.Flank")
maxSum=max(sumBreakdown$Sum)
print('done with IGR')
noIGRHisto <- ggplot(sumBreakdown, aes(x=Variant.Type, y=Sum, fill=Tool, order=-as.numeric(Tool))) + geom_bar(stat="identity", position="dodge") + coord_flip() + geom_text(aes(label=paste(Sum, " of ", Total)), hjust=-0.3, size=3, position = position_dodge(width = 1)) + theme_classic() + guides(fill = guide_legend(reverse = TRUE)) + scale_x_discrete(limits=positions) + scale_y_continuous("Count", expand = c(0, 0), limits=c(0,maxSum+6000))

ggsave(noIGRHisto, filename="HistoNoIGR_4.jpg")



############ Plot FPKM vs sensitivity ###############
library(grid)
library(gtable)
library(gridExtra)

colnames(platSens) <- c("bin", "total", "found", "sensitivity")
colnames(varSens) <- c("bin", "total", "found", "sensitivity")
colnames(mutSens) <- c("bin", "total", "found", "sensitivity")

## Get counts by FPKM bin for not found (sensFPKM) and found (foundsens) variants ##
print(platSens)
platsensFPKM <- as.data.frame(platSens %>% group_by(bin) %>% summarise(totalAll=sum(total)))
platfoundsens <- as.data.frame(platSens %>% group_by(bin) %>% summarise(foundAll=sum(found)))

varsensFPKM <- as.data.frame(varSens %>% group_by(bin) %>% summarise(totalAll=sum(total)))
varfoundsens <- as.data.frame(varSens %>% group_by(bin) %>% summarise(foundAll=sum(found)))

mutsensFPKM <- as.data.frame(mutSens %>% group_by(bin) %>% summarise(totalAll=sum(total)))
mutfoundsens <- as.data.frame(mutSens %>% group_by(bin) %>% summarise(foundAll=sum(found)))

## Join counts of found and not found variants by FPKM bin ##
platsensFPKM <- platsensFPKM %>% left_join(platfoundsens, by="bin")
varsensFPKM <- varsensFPKM %>% left_join(varfoundsens, by="bin")
mutsensFPKM <- mutsensFPKM %>% left_join(mutfoundsens, by="bin")


## Add sensitivity column and replace NAs with 0s ##
platsensFPKM <- as.data.frame(platsensFPKM %>% group_by(bin) %>% mutate(Sensitivity=foundAll / (totalAll)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))
platsensFPKM <- platsensFPKM %>% mutate(Tool="Platypus")
varsensFPKM <- as.data.frame(varsensFPKM %>% group_by(bin) %>% mutate(Sensitivity=foundAll / (totalAll)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))
varsensFPKM <-	varsensFPKM %>% mutate(Tool="VarDict")
mutsensFPKM <- as.data.frame(mutsensFPKM %>% group_by(bin) %>% mutate(Sensitivity=foundAll / (totalAll)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))
mutsensFPKM <- mutsensFPKM %>% mutate(Tool="Mutect2")

## Bind Plat, Var, Mut ##

sensFPKM <- bind_rows(platsensFPKM, varsensFPKM, mutsensFPKM)

print(sensFPKM)

maxCount=max(platsensFPKM$totalAll)
hist_top <- ggplot(platsensFPKM, aes(x=bin, y=totalAll)) + geom_bar(stat="identity") + theme_classic() + scale_x_continuous(limits=c(-13,13), breaks = seq(-13, 13, by = 1)) + scale_y_continuous(trans='log2', breaks=c(5,30,250,2050)) 
#log10(breaks=seq(0,4000, by=1000))
jpeg(filename="sensitivityFPKM_4.jpg", width=10, height=6, units="in", res=300)
hist_top <- hist_top + labs(y="# Known Variants", x="") 
sensPlot <- ggplot(sensFPKM, aes(x=bin, y=Sensitivity, color=Tool, group=Tool)) + geom_point()  + geom_line() + theme_classic() + scale_x_continuous(limits=c(-13,13), breaks = seq(-13, 13, by = 1)) + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by = 0.1))
sensPlot <- sensPlot + labs(x="log2(FPKM)", y="Sensitivity") 

#totalPlot <- grid.arrange(hist_top, sensPlot, nrow=2, heights=c(3,8))
g1 <- ggplotGrob(hist_top)
g2 <- ggplotGrob(sensPlot)
colnames(g1) <- paste0(seq_len(ncol(g1)))
colnames(g2) <- paste0(seq_len(ncol(g2)))

grid.draw(gtable_combine(g1, g2, along=2))
dev.off()
#ggsave(totalPlot, filename="sensitivityFPKM.jpg", width=10, height=6, units="in", dpi=300)




######################## Plot Coverage vs Sensitivity ##################
## concat coverage files for each tool ##
covFound <- bind_rows(platCovFound, varCovFound, mutCovFound)
covNotFound <- bind_rows(platCovNotFound, varCovNotFound, mutCovNotFound)

print(head(covFound))
## set bins by total coverage ## 
# set max read coverage bin to show on plot (goes >15,000 but no variants in most bins) 
binMax <- 450 #max(covNotFound$V7)
print(binMax)
covFound$bin <- cut(covFound$V7, breaks=seq(0, binMax+25, by=25), labels=seq(0, binMax, by=25))
covNotFound$bin <- cut(covNotFound$V7, breaks=seq(0, binMax+25, by=25), labels=seq(0, binMax, by=25))

covFound <- as.data.frame(covFound %>% group_by(Tool, bin) %>% summarise(found=n()))
covNotFound <- as.data.frame(covNotFound %>% group_by(Tool, bin) %>% summarise(notFound=n()))

## join found and not found df's by bin and tool
totCov <- covNotFound %>% full_join(covFound, by=c("bin", "Tool"))
totCov <- totCov %>% drop_na() #mutate_all(funs(ifelse(is.na(.), 0, .)))

## add total counts column
totCov <- as.data.frame(totCov %>% group_by(Tool, bin) %>% mutate(total=as.numeric(found+notFound)))
print("totCov")
print(head(totCov))

## calculate sensitivity and add as column for each bin
totCov <- as.data.frame(totCov %>% group_by(Tool, bin) %>% mutate(sensitivity=(as.numeric(found) / total)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))
print(totCov)

## open jpg file for writing
jpeg(filename="sensitivityCoverage_4.jpg", width=10, height=6, units="in", res=300)
maxCount=max(totCov$total)
histCov <- totCov %>% distinct(totCov, total, .keep_all=TRUE)

## plot histogram of known variants in each bin 
hist_top <- ggplot(histCov, aes(x=bin, y=total)) + geom_bar(stat="identity", position=position_dodge(width=1)) + theme_classic() + scale_y_continuous(trans='log2', breaks=c(5,30,250,2050)) + scale_x_discrete(breaks = seq(0, binMax, by=25)) #+ scale_y_continuous(limits=c(0, maxCount=50), trans='log2')
hist_top <- hist_top + labs(y="# Known Variants", x="")
## plot sensitivity vs read coverage
sensPlot <- ggplot(totCov, aes(x=bin, y=sensitivity, color=Tool, group=Tool)) + geom_point() + geom_line() + theme_classic() + scale_x_discrete(breaks = seq(0, binMax+25, by = 25)) + scale_y_continuous(limits=seq(0,1.0), breaks = seq(0, 1, by = 0.1))
sensPlot <- sensPlot + labs(x="Read Coverage", y="Sensitivity")

## save both plots as grid on same jpg
g1 <- ggplotGrob(hist_top)
g2 <- ggplotGrob(sensPlot)
colnames(g1) <- paste0(seq_len(ncol(g1)))
colnames(g2) <- paste0(seq_len(ncol(g2)))
grid.draw(gtable_combine(g1, g2, along=2))
dev.off()

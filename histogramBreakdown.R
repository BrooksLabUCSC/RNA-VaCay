#!/usr/bin/env Rscript

library("dplyr")
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Supply working directory, whitelist, tumor, normal VCFs, and .n", call.=FALSE)
}


setwd(args[1])
tool=args[2]

breakdown <- read.delim('countBreakdown.txt', sep="\t", header=TRUE, stringsAsFactors=FALSE)
breakdown[is.na(breakdown)] <- 0
averageSensitivity <- as.numeric(breakdown %>% summarize(Mean=mean(Sensitivity)))
averageFiltGermline <- as.numeric(breakdown %>% summarize(Mean=mean(numFiltGermline)))
averageFiltEdits <- as.numeric(breakdown %>% summarize(Mean=mean(numFiltEdits)))
print("Average Sensitivity")
print(averageSensitivity)
print("Average Filtered PON")
print(averageFiltGermline)
print("Average Filtered RNAEdits")
print(averageFiltEdits)

breakdown$Donor.ID <- NULL
breakdown$Sensitivity <- NULL
breakdown$numFiltGermline <- NULL
breakdown$numFiltEdits <- NULL

########## Plot Mean of variants not found (including IGR) ##########

sumBreakdown <- breakdown
breakdown <- breakdown %>% summarize_all(mean)

library(tidyr)

# reshape table
breakdown <- breakdown %>% gather(Variant.Type, Mean) %>% arrange(desc(Mean))

positions <- breakdown$Variant.Type
maxMean <- max(breakdown$Mean)
breakHisto <- ggplot(breakdown, aes(x=Variant.Type, y=Mean)) + geom_bar(stat="identity") + geom_text(aes(label=Mean), hjust=-0.25) + theme_classic() + scale_x_discrete(limits=positions) + scale_y_continuous(expand = c(0, 0), limits=c(0, maxMean+1000)) + coord_flip() + ggtitle(paste("Mean number of variants not found by", tool, sep=" "))

ggsave(breakHisto, filename="meanHisto.jpg")


########## Plot Sum without IGR, 5'Flank, and Introns ###########
sumBreakdown <- sumBreakdown %>% summarize_all(sum)

sumBreakdown$IGR <- NULL
sumBreakdown$Intron <- NULL
sumBreakdown$X5.Flank <- NULL

sumBreakdown <- sumBreakdown %>% gather(Variant.Type, Sum) %>% arrange(desc(Sum))

positions <- sumBreakdown$Variant.Type

maxSum=max(sumBreakdown$Sum)
breakHisto <- ggplot(sumBreakdown, aes(x=Variant.Type, y=Sum)) + geom_bar(stat="identity") + geom_text(aes(label=Sum), hjust=-0.25) + theme_classic() + scale_x_discrete(limits=positions) + scale_y_continuous("Count", expand = c(0, 0), limits=c(0,maxSum+1000)) + coord_flip() + ggtitle(paste("Breakdown of variants not found by", tool, sep=" "))

ggsave(breakHisto, filename="countsHistoNoIGR.jpg")



############ Plot FPKM vs sensitivity ###############
library(gridExtra)

countFPKM <- read.delim('sensitivityFPKM.txt', sep="\t", header=TRUE, stringsAsFactors=FALSE)

colnames(countFPKM) <- c("bin", "total", "found", "sensitivity")
print(countFPKM)
## Get counts by FPKM bin for not found (sensFPKM) and found (foundsens) variants ##
sensFPKM <- as.data.frame(countFPKM %>% group_by(bin) %>% summarise(totalAll=sum(total)))
foundsens <- as.data.frame(countFPKM %>% group_by(bin) %>% summarise(foundAll=sum(found)))

## Join counts of found and not found variants by FPKM bin ##
sensFPKM <- sensFPKM %>% left_join(foundsens, by="bin")

## Add sensitivity column and replace NAs with 0s ##
sensFPKM <- as.data.frame(sensFPKM %>% group_by(bin) %>% mutate(Sensitivity=foundAll / (totalAll)) %>% mutate_all(funs(ifelse(is.na(.), 0, .))))

print(sensFPKM)
maxCount=max(sensFPKM$totalAll)
hist_top <- ggplot(sensFPKM, aes(x=bin, y=totalAll)) + geom_bar(stat="identity") + theme_classic() + scale_x_continuous(limits=c(-13,13), breaks = seq(-13, 13, by = 1)) + scale_y_continuous(trans='log2') 
#log10(breaks=seq(0,4000, by=1000))

hist_top <- hist_top + labs(y="# Known Variants", x="") 
sensPlot <- ggplot(sensFPKM, aes(x=bin, y=Sensitivity)) + geom_point() + theme_classic() + scale_x_continuous(limits=c(-13,13), breaks = seq(-13, 13, by = 1)) + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by = 0.1))
sensPlot <- sensPlot + labs(x="log2(FPKM)", y="Sensitivity") 
totalPlot <- grid.arrange(hist_top, sensPlot, nrow=2, heights=c(3,8))

ggsave(totalPlot, filename="sensitivityFPKM.jpg", width=10, height=6, units="in", dpi=300)
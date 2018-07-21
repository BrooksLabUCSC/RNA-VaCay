library(dplyr)


Maflocs <- whiteMaf %>% semi_join(somaticMinusEdits, by=c("Start_position"="V4"))
mafLocsExp <- DO38901 %>% semi_join(Maflocs, by=c("genes"="Hugo_Symbol"))
filt.DO38 <- DO38901 %>% filter(DO38901, DO38901 > 1)
filtMafExp <- filt.DO38 %>% semi_join(whiteMaf, by=c("genes"="Hugo_Symbol"))
filtMaf <- whiteMaf %>% semi_join(filt.DO38, by=c("Hugo_Symbol"="genes"))

filtMaf %>% group_by(Variant_Classification) %>% summarize(count=n())
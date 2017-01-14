#Install and implement libraries

library(genetics) # load the genetics library
library(dplyr)
library(magrittr)
library(sjPlot)

#Import and attach data
load("EPI556_Data01.rdata")
ls()
attach(fms)
class(fms)

#Using DPLYR to get caucasian observations only
causFMS <- filter(fms, Race == "Caucasian")
finalFMS <- select(causFMS, contains("AKT1"),Gender, Age, NDRM.CH) 
finalFMS[, 1:25] <- sapply(finalFMS[, 1:25], as.numeric)-1

#With Phenotypes

withPheno <- lm(  NDRM.CH~ . , data=finalFMS )


#Without Phenotypes
finalFMS = subset(finalFMS, select = -c(Age,Gender) )

withoutPheno <- lm(  NDRM.CH~ . , data=finalFMS )



sjt.lm(withPheno, withoutPheno)
sjt.lm(withPheno, withoutPheno,
       depvar.labels = c("NDRM.CH With Phenotype", "NDRM.CH Without Phenotype"),
       file = "TableOne")

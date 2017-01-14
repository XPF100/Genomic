load("GenABEL_GWAS1.RData")
ls() 
library(GenABEL)
library(dplyr)

#Clean the data from duplicate rows
SNPNames1=snpnames(GWAS1)
SNPNames2=SNPNames1[!(SNPNames1 %in% c("rs4886982", "rs7965657"))] # remove SNPs from the datasets - 2*2=4 SNP columns 
GWAS2<-GWAS1[phdata(GWAS1)$inferred_population == "Cauc",SNPNames2] 

qc1 <- check.marker(GWAS2, callrate = 0.95, maf = 0.1, perid.call = 0.95, p.level = 0.0001 )

GWAS.qc1 <- GWAS2[qc1$idok,qc1$snpok] # save cleaned samples/SNPs in a new GWAS dataset - GWAS.qc1 


#Fit model and print top 10
model <- qtscore(CYP2D6_activity, GWAS.qc1, trait = "gaussian")
descriptives.scan(model, sort="P1df")[1:10,]


#Calculate Inflation
lambda(model)
estlambda(model[, "P1df"])
estlambda(model[, "P2df"])

### plotting
pval=na.omit(as.numeric(as.character(model[,"P1df"])))
pval.expected=c(1:length(pval))/(length(pval)+1)
plot(-log10(sort(pval.expected)),-log10(sort(pval)),xlim=c(0,5),ylim=c(0,5),cex=0.8,
     xlab="Expected -log(p-value)", ylab="Observed -log(p-value)",pch=19)
abline(a=c(0,0),b=1, col="red")

plot(model)


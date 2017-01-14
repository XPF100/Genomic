mv * /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex



## trait - supar_mean
## For PPG Samples, no people with heart transplant
#Model: supar_mean ~ EnrolledAge + Gender + EPI_CrCl_imputed_binary + PC1-PC10 + snp
vars <- c("supar_rest","supar_45min","supar_90min" )
trait<-c("supar_mean" )
cov<-c("EnrolledAge", "Gender", "EPI_CrCl_imputed_binary") 
setwd("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Chang/PPG_supar_EnrolledAge_Gender_EPI_CrCl_imputed_binary")

#######################################
##### Quantile-Quantile Plot ##########
#######################################
QQplot<-function(p_value,filename,title){
  #Q-Q Plot with CI based on the comments following the article:
  #http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
  #First create the observed and expected P-value vectors
  
  p=-log(na.omit(sort(as.numeric(as.character(p_value)))),10)
  N=length(p)
  p.expected=-log(c(1:N)/(N+1),10)
  MAX=max(c(p, p.expected))
  
  ## create the confidence intervals
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury) or see
  ## http://en.wikipedia.org/wiki/Order_statistic
  
  for(k in 1:N){
    c95[k] <- qbeta(0.95,k,N-k+1)
    c05[k] <- qbeta(0.05,k,N-k+1)
  }
  
  #inflation factor
  
  GC1=dchisq(qchisq(median(sort(na.omit(as.numeric(as.character(p_value))))),df=1),df=1)/dchisq(qchisq(median(c(1:N)/(N+1)),df=1),df=1)
  GC<-format(GC1, digits=4, width=0)
  
  ## plot the two confidence lines and output to .png file
  png(filename = filename, width = 1500, height = 1500,bg = "white",pointsize=24)
  plot(p.expected, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
  par(new=T)
  plot(p.expected, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
  ## add the diagonal
  abline(0,1,col="red")
  par(new=T)
  
  ## add the qqplot
  qqplot(p.expected,p, ylim=c(0,MAX), xlim=c(0,MAX),plot.it=TRUE ,
         xlab="Expected -log(p-value)", ylab="Observed -log(p-value)", pch=19,
         main=paste(title,"\n inflation factor = ",GC,sep="") )
  
  dev.off()
} #End of function "QQplot"

#######################################
########### manhattan plot ############
#######################################
manhattan<-function (result,ID,p,FDR, cbp, pos, main.title,filename)
{
  
  png(filename = filename, width = 1600, height = 900,bg = "white",pointsize=24)
  
  result[,p]=as.numeric(as.character(result[,p]))
  cutoff <- length(result[,p]) - nrow(result[which(result[,FDR]<0.05),]) + 1
  if (cutoff > length(result[,p]))
  { cutoff = length(result[,p])}
  
  
  cbp=result[,cbp]
  pos=result[,pos]
  
  cbp <- as.character(cbp)
  
  cbp[cbp == "X"] = "23"
  cbp[cbp == "Y"] = "24"
  cbp <- as.numeric(cbp)
  
  pos <- as.numeric(pos)
  k = order(cbp, pos)
  cbp = cbp[k]
  pos = pos[k]
  
  score= -log(result[,p][k], base = 10)
  cmax = max(cbp)
  genomepos = pos
  
  
  
  cbp_unique <- unique(cbp)
  
  #give a 100 interval between each cbp
  if (length(cbp_unique) > 1) {
    for (r in cbp_unique[2:length(cbp_unique)]) {
      genomepos[cbp == r] = genomepos[cbp == r] + max(genomepos[cbp == cbp_unique[which(cbp_unique == r) - 1]]) + 100
    }
  }
  bonsig <- which(score > -log(0.05/cutoff, base = 10))
  lessig <- which(score < -log(0.05/cutoff, base = 10) & score > -log(0.05/cutoff, base = 10)/2)
  
  cbpcolor <- cbp
  original <- par(no.readonly = TRUE)
  
  
  
  cbpcolor[which(cbpcolor == "7" | cbpcolor == "15" | cbpcolor == "23")] = "orange"
  
  
  
  y.lab <- c(0, max(-log(0.05/cutoff, base = 10), score, na.rm = TRUE) +  0.05)
  plot(range(genomepos), y.lab, type = "n", xaxt = "n", bty = "7",
       xlab = "cbpomosome", ylab = expression(paste("Observed -log ",scriptstyle(10), "(P-values)", sep = "")), main = main.title)
  points(genomepos, score, col = cbpcolor, cex = 0.2)
  points(genomepos[bonsig], score[bonsig], col = cbpcolor[bonsig], pch = 16)
  points(genomepos[lessig], score[lessig], col = cbpcolor[lessig], pch = 16, cex = 0.5)
  abline(-log(0.05/cutoff, base = 10), 0,col="red")
  abline(5,0,col="orange" )
  
  
  # if there is a significant FDR
  temp= as.numeric(as.character(result[,FDR]))
  
  if (nrow(result[which(temp <0.05),]) > 0) {
    #draw a line which is the largest significant p
    abline(-log(max(as.numeric(as.character(result[which(temp<0.05),p]))), base = 10), 0,col="green", lty = 2)
    
    if(FALSE){
      FDR.cutoff <- -log(max(as.numeric(as.character(result[which(temp<0.05),p]))), base = 10)
      index=result[,ID]
      FDR.labels <- index[k][which(score > FDR.cutoff)]
      text(genomepos[which(score > FDR.cutoff)], score[which(score > FDR.cutoff)], labels = FDR.labels, pos = "4",cex = 0.7)
    }
  }
  
  
  meds = matrix(0, cmax)
  labs = matrix(NA, cmax)
  tix = matrix(0, cmax + 1)
  labs_null = matrix(NA, cmax + 1)
  
  for (t in 1:length(cbp_unique)) {
    tix[t + 1] = max(genomepos[cbp == cbp_unique[t]])
    meds[t] <- (tix[t] + tix[t + 1, 1])/2
    labs[t] = cbp_unique[t]
  }
  
  if (cmax >= 23) {
    labs[which(labs == 23)] = "X"
    if (sum(cbp == 24) > 0) {
      labs[which(labs == 24)] = "Y"
    }
  }
  axis(side = 1, at = meds, lwd.ticks = 0, labels = labs)
  axis(side = 1, at = tix, lwd.ticks = 1, labels = labs_null)
  
  dev.off()
  
}  #End of function "manhattan"
## Definition of functions finished

#####################
####Calculate MAF####
#####################
#Open another window and use plink. I copied the PLINK files into my directory so that we keep output files in it
#export PATH=$PATH:/home/qhui/bin
#cd /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Chang/PPG_supar_EnrolledAge_Gender_EPI_CrCl_imputed_binary/Alex
#plink --bfile SNP_Sample_Filtered_MIPS_MIMS_IdenticalIDFixed_NoNa_NoDup_NoRaceMisMatch_AA --freq --out PPG_MAF_AA --noweb
#plink --bfile SNP_Sample_Filtered_MIPS_MIMS_IdenticalIDFixed_NoNa_NoDup_NoRaceMisMatch_WT --freq --out AAPPG_MAF_WT --noweb

#####################################Code for trait interested starts here####################################
#Model: supar_mean ~ EnrolledAge + Gender + EPI_CrCl_imputed_binary + PC1-PC10 + snp

setwd("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Chang/PPG_supar_EnrolledAge_Gender_EPI_CrCl_imputed_binary/")

phenos1<-get(load(file="/sun_group/sunlab1/projects/PPG_GWAS/Data/MIPS_suPAR.rdata"))
dim(phenos1)  #699 198
names(phenos1)
#write.csv(phenos1, file=paste("phenos1.csv", sep=""), row.names=F, quote=F)
#"AfricanAmerican" "Caucasian"
#"supar_rest" "supar_45min" "supar_90min"
#"ID"(1001, 1002) "EnrolledAge" "Gender"(1-Female, 2-Male) "SmokeHistory" "BMI" "CAD_MI" "HeartFailure" "Diabetes" "Hypertension" "Dyslipidemia" "EPI_CrCl_imputed" 
phenos1<-phenos1[,c("ID","supar_rest","supar_45min","supar_90min","EnrolledAge","Gender","SmokeHistory","BMI","CAD_MI","HeartFailure","Diabetes","Hypertension","Dyslipidemia","EPI_CrCl_imputed")]

phenos2<-get(load(file="/sun_group/sunlab1/projects/PPG_GWAS/Data/MIPS_9_19.rdata"))
dim(phenos2)   #695 1710
names(phenos2)
#write.csv(phenos2, file=paste("phenos2.csv", sep=""), row.names=F, quote=F)
#"AfricanAmerican" "Caucasian"
#"ID"(1001, 1002) "EnrolledAge" "Gender"(1-Female, 2-Male) "SmokeHistory" "BMI" "CAD_MI" "HeartFailure" "Diabetes" "Hypertension" "Dyslipidemia" "EPI_CrCl_imputed" "CAD50binary"
phenos2<-phenos2[,c("ID","AfricanAmerican","Caucasian","CAD50binary")]

phenos<-merge(phenos1, phenos2, by="ID")
# 698  17
names(phenos)
#[1] "ID"               "supar_rest"       "supar_45min"      "supar_90min"     
#[5] "EnrolledAge"      "Gender"           "SmokeHistory"     "BMI"             
#[9] "CAD_MI"           "HeartFailure"     "Diabetes"         "Hypertension"    
#[13] "Dyslipidemia"     "EPI_CrCl_imputed" "AfricanAmerican"  "Caucasian"       
#[17] "CAD50binary"     
names(phenos)[names(phenos)=="ID"] <- "IID"
phenos$IID<-paste("MIPS", phenos$IID, sep="")

phenos$supar_mean <- rowMeans(phenos[,vars], na.rm=TRUE)
phenos <- subset(phenos, phenos$supar_mean >= 0)
phenos$supar_mean <- log(phenos$supar_mean)
names(phenos)

###################################GOOD TO HERE################################


####################For WT#######################

PCs_WT<-read.csv(file="/sun_group/sunlab1/projects/PPG_GWAS/Work_11032015/pca_10_PPG_NoRaceMisMatch_WT.csv", header=T, as.is=T)
dim(PCs_WT) #375 11
#"IID"  "PC1"  "PC2"  "PC3"  "PC4"  "PC5"  "PC6"  "PC7"  "PC8"  "PC9"  "PC10"

#Read the MAF for WT
MAF_WT<-read.table(file="/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/AAPPG_MAF_WT.frq", sep="", header=T, as.is=T)
dim(MAF_WT)     #1614035       6
#MAF_WT[1:3,]
#    CHR             SNP A1 A2     MAF NCHROBS
#1   0 1:109819296-C-A  0  A 0.00000     750
#2   0 1:110229944-G-C  G  C 0.04545     726
#3   0  1:11082347-G-C  0  G 0.00000     750

phenos_WT<-phenos[phenos$IID %in% substring(PCs_WT$IID, 1, 8),]  #357 17
phenos_WT$IID<-ifelse(phenos_WT$IID %in% setdiff(phenos_WT$IID, PCs_WT$IID), paste(phenos_WT$IID, "-D00", sep=""), phenos_WT$IID)
setdiff(phenos_WT$IID, PCs_WT$IID)
setdiff(PCs_WT$IID, phenos_WT$IID) #All MIMSxxxx
#Only kept the MIPSxxxx that were matched in PCs_WT. The MIMSxxxx in PCs_WT were excluded.
dim(phenos_WT) #357 17

#install.packages("SKAT")
library(SKAT)
FamilyID<-Read_Plink_FAM("/sun_group/sunlab1/projects/PPG_GWAS/Work_11032015/SNP_Sample_Filtered_MIPS_MIMS_IdenticalIDFixed_NoNa_NoDup_NoRaceMisMatch_WT.fam")
dim(FamilyID)    #375   6
FamilyID<-FamilyID[,c("FID","IID")]

#merge the phenos, PCs and FID by IID
phenos_PCs_WT <- merge(phenos_WT,PCs_WT, by="IID")
dim(phenos_PCs_WT)
#357 27
phenos_PCs_WT <- merge(FamilyID,phenos_PCs_WT, by="IID")
phenos_PCs_WT$EPI_CrCl_imputed_binary<-ifelse(phenos_PCs_WT$EPI_CrCl_imputed>=60, 1, 0)
#NA remains NA
dim(phenos_PCs_WT)   
#357 29
names(phenos_PCs_WT)
#[1] "IID"              "FID"              "supar_rest"       "supar_45min"     
#[5] "supar_90min"      "EnrolledAge"      "Gender"           "SmokeHistory"    
#[9] "BMI"              "CAD_MI"           "HeartFailure"     "Diabetes"        
#[13] "Hypertension"     "Dyslipidemia"     "EPI_CrCl_imputed" "AfricanAmerican" 
#[17] "Caucasian"        "CAD50binary"      "PC1"              "PC2"             
#[21] "PC3"              "PC4"              "PC5"              "PC6"             
#[25] "PC7"              "PC8"              "PC9"              "PC10"           "EPI_CrCl_imputed_binary" 

#Check distribution of the dependent variable.
######Before turning NAs into "-9"!!!!
dim(phenos_PCs_WT) # 357 29
summary(phenos_PCs_WT$supar_mean) #117 missing
png(filename =paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_Hist_supar_mean_WT.png", sep=""), width = 1500, height = 1500,bg = "white",pointsize=24)
hist(phenos_PCs_WT$supar_mean, col="green")
dev.off()

summary(phenos_PCs_WT$EnrolledAge) #1 missing
png(filename =paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_Hist_EnrolledAge_WT.png", sep=""), width = 1500, height = 1500,bg = "white",pointsize=24)
hist(phenos_PCs_WT$EnrolledAge, col="green")
dev.off()

summary(phenos_PCs_WT$Gender) #1 missing
table(phenos_PCs_WT$Gender)
#1  2 
#80 276
#(1-Female, 2-Male)

summary(phenos_PCs_WT$EPI_CrCl_imputed) #1 missing
png(filename =paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_Hist_EPI_CrCl_imputed_WT.png", sep=""), width = 1500, height = 1500,bg = "white",pointsize=24)
hist(phenos_PCs_WT$EPI_CrCl_imputed, col="green")
dev.off()

summary(phenos_PCs_WT$EPI_CrCl_imputed_binary) #1 missing
table(phenos_PCs_WT$EPI_CrCl_imputed_binary)
#0   1 
#64 292 

outcome_WT <- phenos_PCs_WT[,c("FID", "IID","supar_mean")] 
dim(outcome_WT)    #357  3
outcome_WT[is.na(outcome_WT)] <- "-9"   
#if assign the value as numeric -9, it's gonna be -9.000, which may not be recognized by plink

#pick the covars
cov_WT<-phenos_PCs_WT[, c("FID","IID", cov, "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
dim(cov_WT)
#357  15 
cov_WT[is.na(cov_WT)] <- "-9"

n=1

#output the two files for PLINK
write.table(outcome_WT, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_outcome_WT_EPI_CrCl_imputed_binary_", trait[n], ".txt", sep=""), sep=" ", row.names=F, quote=F)
write.table(cov_WT, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_cov_WT_EPI_CrCl_imputed_binary_", trait[n], ".txt", sep=""), sep=" ", row.names=F, quote=F)

########################################################
####### Run association calculation using PLINK ########
########################################################
##Method 1, use command line commands: (suggestion: open another Putty window and use screen)
export PATH=$PATH:/home/qhui/bin
cd /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Chang/PPG_supar_EnrolledAge_Gender_EPI_CrCl_imputed_binary
plink --bfile SNP_Sample_Filtered_MIPS_MIMS_IdenticalIDFixed_NoNa_NoDup_NoRaceMisMatch_WT --pheno /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_outcome_WT_EPI_CrCl_imputed_binary_supar_mean.txt --pheno-name supar_mean --covar /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_cov_WT_EPI_CrCl_imputed_binary_supar_mean.txt --covar-name EnrolledAge,Gender,EPI_CrCl_imputed_binary,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --linear --out /sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/PPG_WT_EPI_CrCl_imputed_binary_supar_mean --noweb
######################################

###########Good#######################

#If you used method 1, go back to the Putty window where you are running R1
setwd("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex")
results_supar_mean_WT<-read.table(file=paste("PPG_WT_EPI_CrCl_imputed_binary_supar_mean.assoc.linear", sep=""), sep="", header=T, as.is=T)
#This can take a while

dim(results_supar_mean_WT)
#22641840        9
#Each snp has 14 rows
#But the number of rows is not a multiply of 14 for some reason.
#Removing non-autosomal chromosoms seems resolved this puzzle:
results_supar_mean_WT_Auto<- results_supar_mean_WT[results_supar_mean_WT$CHR>0 & results_supar_mean_WT$CHR<23,]
dim(results_supar_mean_WT_Auto)
#21609896        9
#21609896/14 = 1543564

#Save the results as R data so it can be accessed quicker when needed
save(results_supar_mean_WT_Auto, file=paste("PPG_WT_EPI_CrCl_imputed_binary_supar_mean_Auto.assoc.linear.rda", sep=""),compress=TRUE, compression_level=9)

saveRDS(results_supar_mean_WT_Auto, file=paste("PPG_WT_EPI_CrCl_imputed_binary_supar_mean_Auto.assoc.linear.rda", sep=""))

#This is how the first a few rows look like:
results_supar_mean_WT_Auto[1:4,]
#      CHR         SNP    BP A1                    TEST NMISS BETA STAT  P
#178095   1 JHU_1.17537 17538  0                     ADD   238   NA   NA NA
#178096   1 JHU_1.17537 17538  0             EnrolledAge   238   NA   NA NA
#178097   1 JHU_1.17537 17538  0                  Gender   238   NA   NA NA
#178098   1 JHU_1.17537 17538  0 EPI_CrCl_imputed_binary   238   NA   NA NA

#The row with TEST=ADD is the additive effect, which is what we need
WT_ADD<-results_supar_mean_WT_Auto[results_supar_mean_WT_Auto$TEST=="ADD",]
dim(WT_ADD)
# 1543564        9

#Calculate the multi-comparison-corrected p-values and add to the data frame
WT_ADD$BH=p.adjust(WT_ADD[,"P"], method = "BH", n = length(WT_ADD[,"P"]))  ## FDR Adjustment. The n in this line is not the n for trait loop but the total number of SNPs
WT_ADD$bonferroni=p.adjust(WT_ADD[,"P"], method = "bonferroni", n = length(WT_ADD[,"P"]))  ##Bonferroni Ajustment

#Prepare file for META analysis
#Get MAF and allel names for each SNP
WT_ADD_AUTO_MAF<-merge(WT_ADD, MAF_WT[, c("SNP", "MAF", "A2")], by="SNP", all.x=T, sort=F)  #WT_ADD already has 1st allel, A1

#As PLINK output does not include SE, which is needed for META analysis, we calculate it
WT_ADD_AUTO_MAF$SE<-ifelse(is.na(WT_ADD_AUTO_MAF$BETA), NA, WT_ADD_AUTO_MAF$BETA/WT_ADD_AUTO_MAF$STAT)

WT_ADD_AUTO_MAF<-WT_ADD_AUTO_MAF[, c("CHR", "SNP", "BP", "A1", "A2", "MAF", "NMISS", "BETA", "SE", "STAT", "P", "BH", "bonferroni")]
dim(WT_ADD_AUTO_MAF) 
# 1543564      13

n=1
write.table(WT_ADD_AUTO_MAF, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_WT_EPI_CrCl_imputed_binary_ADD_AUTO_MAF_", trait[n],".txt", sep=""), sep=" ", row.names=F, quote=F)   #This is the file for META Analysis input

#Prepare data for QQ plot
p_value<-WT_ADD_AUTO_MAF[!is.na(WT_ADD_AUTO_MAF$P), "P"]

filename=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_QQ_WT_EPI_CrCl_imputed_binary_", trait[n], ".png", sep="")
title=paste("Q-Q Plot of PPG WT - EnrolledAge, Gender, EPI_CrCl_imputed_binary, ", trait[n], sep="")

#Call the function defined above to draw QQ plot
QQplot(p_value,filename,title)

#The QQ plot shows a "stage" shape
#It's caused by rear snps. 
#Remove the SNPs with MAF lower than 5% and try again:
WT_ADD_AUTO_COMMON<-WT_ADD_AUTO_MAF[WT_ADD_AUTO_MAF$MAF>0.05,]
dim(WT_ADD_AUTO_COMMON)
#  605531     13

p_value<-WT_ADD_AUTO_COMMON[!is.na(WT_ADD_AUTO_COMMON$P), "P"]

filename=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_QQ_WT_EPI_CrCl_imputed_binary_", trait[n], "_AUTO_MAF05.png", sep="")
title=paste("Q-Q Plot of PPG WT - EnrolledAge, Gender, EPI_CrCl_imputed_binary, ", trait[n], ", Auto SNP MAF>0.05", sep="")

#Call the function defined above to draw QQ plot
QQplot(p_value,filename,title)
#This time QQ looks good

# Prepare data for Manhattan plot
WT_ADD_AUTO_COMMON$Outcome<-paste(trait[n], sep="")

ID="Outcome"
p="P"
FDR="BH"
cbp="CHR"
pos="BP"
main.title=paste("EnrolledAge, Gender, EPI_CrCl_imputed_binary - \nManhattan plot of the association between PPG WT SNPs and ", trait[n], sep="")
filename=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_Manhattan_PPG_WT_EPI_CrCl_imputed_binary_", trait[n], "_MAF.png", sep="")

#Call the function defined above to draw Manhattan plot
manhattan(WT_ADD_AUTO_COMMON,ID,p,FDR, cbp, pos, main.title, filename)

#Make saummary table
#First a simple table containing SNPs with p-value smaller than 1e-5
WT_Sig_COMMON<-WT_ADD_AUTO_COMMON[!is.na(WT_ADD_AUTO_COMMON$P) & WT_ADD_AUTO_COMMON$P<1e-5, ]
dim(WT_Sig_COMMON)
# 185 14
write.csv(WT_Sig_COMMON, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_WT_EPI_CrCl_imputed_binary_", trait[n], "_Sig_MAF05.csv", sep=""), row.names=F, quote=F)

#To include covariates in the summary table, except PCs
WT_Sig_COMMON_Cov<-results_supar_mean_WT_Auto[results_supar_mean_WT_Auto$SNP %in% WT_Sig_COMMON$SNP & results_supar_mean_WT_Auto$TEST %in% c(cov), ]
dim(WT_Sig_COMMON_Cov)   
# 555    9

WT_Sig_COMMON_Complete<-matrix(0, dim(WT_Sig_COMMON)[1], (9+3*(length(cov)+1)))
dim(WT_Sig_COMMON_Complete)
# 185  21

WT_Sig_COMMON_Complete<-as.data.frame(WT_Sig_COMMON_Complete)
colnames(WT_Sig_COMMON_Complete)[1:7]<-c("Outcome","SNP","CHR","BP","A1","NMISS","MAF")
colnames(WT_Sig_COMMON_Complete)[(9+3*(length(cov)+1)-4):(9+3*(length(cov)+1))]<-c("BETA_SNP", "STAT_SNP","P_SNP", "BH", "bonferroni")

#The following loop will assign values for each cell, and also put the column names for covariates.
for(i in 1:dim(WT_Sig_COMMON)[1]){
  WT_Sig_COMMON_Complete[i,1:7]<-WT_Sig_COMMON[i,c("Outcome","SNP","CHR","BP","A1","NMISS","MAF")]
  
  for(j in 1:(length(cov))){
    WT_Sig_COMMON_Complete[i,(8+(j-1)*3):(10+(j-1)*3)]<-WT_Sig_COMMON_Cov[(i-1)*3+j, c("BETA","STAT","P")]   
    
    colnames(WT_Sig_COMMON_Complete)[(8+(j-1)*3):(10+(j-1)*3)]<-c(paste("BETA_",cov[j],sep=""), paste("STAT_",cov[j],sep=""),paste("P_",cov[j],sep=""))
  } #end of j loop      
  WT_Sig_COMMON_Complete[i,(9+3*(length(cov)+1)-4):(9+3*(length(cov)+1))]<-WT_Sig_COMMON[i, c("BETA","STAT","P",  "BH", "bonferroni")]
}   #end of i loop

#re-order the rows so the SNP with smallest p-value is on the top.
WT_Sig_COMMON_Complete<-WT_Sig_COMMON_Complete[order(as.numeric(as.character(WT_Sig_COMMON_Complete$P_SNP))),]

write.csv(WT_Sig_COMMON_Complete, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_Summary_PPG_WT_EPI_CrCl_imputed_binary_", WT_Sig_COMMON_Complete$Outcome[1], "_MAF_WithCov.csv", sep=""), row.names=F, quote=F)

#To draw regional plot for the most significant snp, output the corresponding chromosome
#I couldn't do this becasue the most significant SNP does not have a reference SNP name that can be regognized
reg<-WT_ADD_AUTO_COMMON[WT_ADD_AUTO_COMMON$CHR==WT_Sig_COMMON_Complete$CHR[1], c("SNP", "P")]
colnames(reg)<-c("MarkerName", "P.value")
write.table(reg, file=paste("/sun_group/sunlab1/projects/ECRI_Biobank_GWAS/Work/Alex/Alex_PPG_WT_EPI_CrCl_imputed_binary_", WT_Sig_COMMON_Complete$Outcome[1], "_CHR", WT_Sig_COMMON_Complete$CHR[1], "_SNP_", WT_Sig_COMMON_Complete$SNP[1],".txt", sep=""), sep="\t", row.names=F, quote=F)
#Upload this list to http://locuszoom.sph.umich.edu/locuszoom/genform.php?type=yourdata to draw regional plot


#######################################META##########################################
#Then run the Meta analysis using "metal", for black and white, respectively. Example commands can be found in PPG_META_summary_bp.r
#Then use PPG_META_summary_bp.r to summarize
11111111

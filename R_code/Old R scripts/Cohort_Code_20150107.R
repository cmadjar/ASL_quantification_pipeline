# clear workspace
rm(list = ls())
library(ggplot2)
library(labeling)
library(reshape2)


# load datasets
data2	= read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/QC/for analyses/ALL_DATA_20150105.csv",h=T,sep=",")
# convert factorial data to binary data
data2$DCCID						 = as.factor(data2$DCCID)
data2$FH_code					 = as.factor(data2$FH_code)   
data2$Gender_bin					 = as.factor(data2$Gender_bin)
data2$ApoE   					 = as.factor(data2$ApoE)
data2$Apoe_cat  					 = as.factor(data2$Apoe_cat)
data2$K_variant_positive 		 = as.factor(data2$K_variant_positive)
data2$BDNF_allele_no				 = as.factor(data2$BDNF_allele_no)
data2$BDNF_copie_no				 = as.factor(data2$BDNF_copie_no)
data2$Intron_M_allele_no 		 = as.factor(data2$Intron_M_allele_no)
data2$Intron_M_Protective_Variant = as.factor(data2$Intron_M_Protective_Variant) 
data2$ptau_BL_bin_median			 = as.factor(data2$ptau_BL_bin_median)
data2$E4_allele_no				 = as.factor(data2$E4_allele_no)
 
# display the headers and a summary
names(data2)
summary(data2)


# create table with up till 24 month cohort pass CBF data only
cohort_FU24_pass = subset(data2, Project=="PRE" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU12=="pass" & Precuneus_QC_FU12=="pass" & PCC_QC_FU12=="pass" & ACC_QC_FU12=="pass" & GM_QC_FU24=="pass" & Precuneus_QC_FU24=="pass" & PCC_QC_FU24=="pass" & ACC_QC_FU24=="pass")
# remove non used column
cohort_FU24_pass$ACC_Avg_FU03 <- NULL
cohort_FU24_pass$ACC_QC_FU03 <- NULL
cohort_FU24_pass$PCC_Avg_FU03 <- NULL
cohort_FU24_pass$PCC_QC_FU03 <- NULL
cohort_FU24_pass$Precuneus_Avg_FU03 <- NULL
cohort_FU24_pass$Precuneus_QC_FU03 <- NULL
cohort_FU24_pass$GM_Avg_FU03 <- NULL
cohort_FU24_pass$GM_QC_FU03 <- NULL
cohort_FU24_pass_long = cohort_FU24_pass[,c(1, 7, 10, 15, 17, 20, 23, 31, 32, 33, 34, 35, 36, 37, 38, 39, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83)]

#### GM CBF
## organize data (grep GM CBF only and convert table to a longitudinal one)
cohort_FU24_pass_GM_CBF = cohort_FU24_pass_long[,c(1,2,3,4,5,6,7,17,21,25)]
GM_CBF_long = melt(cohort_FU24_pass_GM_CBF, id=c(1,2,3,4,5,6,7), value.name="GM_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_GM <- ggplot(data=GM_CBF_long, aes(x= TimePoint, y= GM_Avg_CBF, group=DCCID))
# spaghetti plot with mean only (could be changed to median if wanted to)
# p + geom_line() + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)
# spaghetti plot with quantile
# p + geom_line() + stat_summary(aes(group=1), geom="point", fun.y=quantile, probs=c(0.25,0.75), shape = 17, size =3) + facet_grid(. ~ E4_allele_no)
# spaghetti plot with mean and standard error shading
p_GM + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM Precuneus
## organize data
cohort_FU24_pass_Prec_CBF = cohort_FU24_pass_long[,c(1,2,3,4,5,6,7,18,22,26)]
Prec_CBF_long = melt(cohort_FU24_pass_Prec_CBF, id=c(1,2,3,4,5,6,7), value.name="Prec_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_Prec <- ggplot(data=Prec_CBF_long, aes(x= TimePoint, y= Prec_Avg_CBF, group=DCCID))
p_Prec + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM PCC
## organize data
cohort_FU24_pass_PCC_CBF = cohort_FU24_pass_long[,c(1,2,3,4,5,6,7,19,23,27)]
PCC_CBF_long = melt(cohort_FU24_pass_PCC_CBF, id=c(1,2,3,4,5,6,7), value.name="PCC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_PCC <- ggplot(data=PCC_CBF_long, aes(x= TimePoint, y=PCC_Avg_CBF, group=DCCID))
p_PCC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM ACC
## organize data
cohort_FU24_pass_ACC_CBF = cohort_FU24_pass_long[,c(1,2,3,4,5,6,7,20,24,28)]
ACC_CBF_long = melt(cohort_FU24_pass_ACC_CBF, id=c(1,2,3,4,5,6,7), value.name="ACC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_ACC <- ggplot(data=ACC_CBF_long, aes(x= TimePoint, y=ACC_Avg_CBF, group=DCCID))
p_ACC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)





######## create table with 2 visits pass CBF data only
# create table with up till 12 month cohort pass CBF data only
cohort_FU12_pass = subset(data2, Project=="PRE" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU12=="pass" & Precuneus_QC_FU12=="pass" & PCC_QC_FU12=="pass" & ACC_QC_FU12=="pass")
# remove non used column
cohort_FU12_pass$ACC_Avg_FU03 <- NULL
cohort_FU12_pass$ACC_QC_FU03 <- NULL
cohort_FU12_pass$PCC_Avg_FU03 <- NULL
cohort_FU12_pass$PCC_QC_FU03 <- NULL
cohort_FU12_pass$Precuneus_Avg_FU03 <- NULL
cohort_FU12_pass$Precuneus_QC_FU03 <- NULL
cohort_FU12_pass$GM_Avg_FU03 <- NULL
cohort_FU12_pass$GM_QC_FU03 <- NULL
cohort_FU12_pass$ACC_Avg_FU24 <- NULL
cohort_FU12_pass$ACC_QC_FU24 <- NULL
cohort_FU12_pass$PCC_Avg_FU24 <- NULL
cohort_FU12_pass$PCC_QC_FU24 <- NULL
cohort_FU12_pass$Precuneus_Avg_FU24 <- NULL
cohort_FU12_pass$Precuneus_QC_FU24 <- NULL
cohort_FU12_pass$GM_Avg_FU24 <- NULL
cohort_FU12_pass$GM_QC_FU24 <- NULL
cohort_FU12_pass_long = cohort_FU12_pass[,c(1, 7, 10, 15, 17, 20, 23, 42, 44, 46, 48, 50, 52, 54, 56, 58, 59, 60, 61, 62, 63, 70, 71, 72, 73, 74, 75)]

#### GM CBF
## organize data (grep GM CBF only and convert table to a longitudinal one)
cohort_FU12_pass_GM_CBF = cohort_FU12_pass_long[,c(1,2,3,4,5,6,8,12)]
GM_CBF_long = melt(cohort_FU12_pass_GM_CBF, id=c(1,2,3,4,5,6), value.name="GM_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_GM <- ggplot(data=GM_CBF_long, aes(x= TimePoint, y= GM_Avg_CBF, group=DCCID))
p_GM + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM Precuneus 
cohort_FU12_pass_Prec_CBF = cohort_FU12_pass_long[,c(1,2,3,4,5,6,9,13)]
Prec_CBF_long = melt(cohort_FU12_pass_Prec_CBF, id=c(1,2,3,4,5,6), value.name="Prec_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_Prec <- ggplot(data=Prec_CBF_long, aes(x= TimePoint, y=Prec_Avg_CBF, group=DCCID))
p_Prec + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM PCC 
cohort_FU12_pass_PCC_CBF = cohort_FU12_pass_long[,c(1,2,3,4,5,6,10,14)]
PCC_CBF_long = melt(cohort_FU12_pass_PCC_CBF, id=c(1,2,3,4,5,6), value.name="PCC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_PCC <- ggplot(data=PCC_CBF_long, aes(x= TimePoint, y=PCC_Avg_CBF, group=DCCID))
p_PCC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM ACC
cohort_FU12_pass_ACC_CBF = cohort_FU12_pass_long[,c(1,2,3,4,5,6,11,15)]
ACC_CBF_long = melt(cohort_FU12_pass_ACC_CBF, id=c(1,2,3,4,5,6), value.name="ACC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_ACC <- ggplot(data=ACC_CBF_long, aes(x= TimePoint, y=ACC_Avg_CBF, group=DCCID))
p_ACC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)




# clear workspace
rm(list = ls())
library(ggplot2)
library(labeling)

# load datasets
data	= read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/QC/for analyses/ALL_DATA_20150105.csv",h=T,sep=",")
# convert factorial data to binary data
data$DCCID						 = as.factor(data$DCCID)
data$FH_code					 = as.factor(data$FH_code)   
data$Gender_bin					 = as.factor(data$Gender_bin)
data$ApoE   					 = as.factor(data$ApoE)
data$Apoe_cat  					 = as.factor(data$Apoe_cat)
data$K_variant_positive 		 = as.factor(data$K_variant_positive)
data$BDNF_allele_no				 = as.factor(data$BDNF_allele_no)
data$BDNF_copie_no				 = as.factor(data$BDNF_copie_no)
data$Intron_M_allele_no 		 = as.factor(data$Intron_M_allele_no)
data$Intron_M_Protective_Variant = as.factor(data$Intron_M_Protective_Variant) 
data$ptau_BL_bin_median			 = as.factor(data$ptau_BL_bin_median)
data$E4_allele_no				 = as.factor(data$E4_allele_no)
 
# display the headers and a summary
names(data)
summary(data)


# create tables with 3 LPs NAP pass CBF data only
nap_3LP_pass	 = subset(data, Project=="NAP" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU03=="pass" & Precuneus_QC_FU03=="pass" & PCC_QC_FU03=="pass" & ACC_QC_FU03=="pass" & GM_QC_FU12=="pass" & Precuneus_QC_FU12=="pass" & PCC_QC_FU12=="pass" & ACC_QC_FU12=="pass")
# remove non used column
nap_3LP_pass$ACC_Avg_FU24 <- NULL
nap_3LP_pass$ACC_QC_FU24 <- NULL
nap_3LP_pass$PCC_Avg_FU24 <- NULL
nap_3LP_pass$PCC_QC_FU24 <- NULL
nap_3LP_pass$Precuneus_Avg_FU24 <- NULL
nap_3LP_pass$Precuneus_QC_FU24 <- NULL
nap_3LP_pass$GM_Avg_FU24 <- NULL
nap_3LP_pass$GM_QC_FU24 <- NULL
nap_3LP_pass_long = nap_3LP_pass[,c(1, 7, 10, 15, 17, 20, 23, 31, 32, 33, 34, 35, 36, 37, 38, 39, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83)]

#### GM CBF
## organize data (grep GM CBF only and convert table to a longitudinal one)
NAP_3LP_pass_GM_CBF = nap_3LP_pass_long[,c(1,2,3,4,5,6,7,17,21,25)]
GM_CBF_long = melt(NAP_3LP_pass_GM_CBF, id=c(1,2,3,4,5,6,7), value.name="GM_Avg_CBF", variable.name="TimePoint")
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
NAP_3LP_pass_Prec_CBF = nap_3LP_pass_long[,c(1,2,3,4,5,6,7,18,22,26)]
Prec_CBF_long = melt(NAP_3LP_pass_Prec_CBF, id=c(1,2,3,4,5,6,7), value.name="Prec_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_Prec <- ggplot(data=Prec_CBF_long, aes(x= TimePoint, y= Prec_Avg_CBF, group=DCCID))
p_Prec + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM PCC
## organize data
NAP_3LP_pass_PCC_CBF = nap_3LP_pass_long[,c(1,2,3,4,5,6,7,19,23,27)]
PCC_CBF_long = melt(NAP_3LP_pass_PCC_CBF, id=c(1,2,3,4,5,6,7), value.name="PCC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_PCC <- ggplot(data=PCC_CBF_long, aes(x= TimePoint, y=PCC_Avg_CBF, group=DCCID))
p_PCC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM ACC
## organize data
NAP_3LP_pass_ACC_CBF = nap_3LP_pass_long[,c(1,2,3,4,5,6,7,20,24,28)]
ACC_CBF_long = melt(NAP_3LP_pass_ACC_CBF, id=c(1,2,3,4,5,6,7), value.name="ACC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_ACC <- ggplot(data=ACC_CBF_long, aes(x= TimePoint, y=ACC_Avg_CBF, group=DCCID))
p_ACC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)





######## create table with 2 LPs pass CBF data only
nap_2LP_pass	 = subset(data, Project=="NAP" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU03=="pass" & Precuneus_QC_FU03=="pass" & PCC_QC_FU03=="pass" & ACC_QC_FU03=="pass")
# remove non used column
nap_2LP_pass$ACC_Avg_FU24 <- NULL
nap_2LP_pass$ACC_QC_FU24 <- NULL
nap_2LP_pass$PCC_Avg_FU24 <- NULL
nap_2LP_pass$PCC_QC_FU24 <- NULL
nap_2LP_pass$Precuneus_Avg_FU24 <- NULL
nap_2LP_pass$Precuneus_QC_FU24 <- NULL
nap_2LP_pass$GM_Avg_FU24 <- NULL
nap_2LP_pass$GM_QC_FU24 <- NULL
nap_2LP_pass$ACC_Avg_FU12 <- NULL
nap_2LP_pass$ACC_QC_FU12 <- NULL
nap_2LP_pass$PCC_Avg_FU12 <- NULL
nap_2LP_pass$PCC_QC_FU12 <- NULL
nap_2LP_pass$Precuneus_Avg_FU12 <- NULL
nap_2LP_pass$Precuneus_QC_FU12 <- NULL
nap_2LP_pass$GM_Avg_FU12 <- NULL
nap_2LP_pass$GM_QC_FU12 <- NULL
nap_2LP_pass$RBANS_immed_mem_FU12 <- NULL
nap_2LP_pass$RBANS_visuo_const_FU12 <- NULL
nap_2LP_pass$RBANS_delay_mem_FU12 <- NULL
nap_2LP_pass$RBANS_language_FU12 <- NULL
nap_2LP_pass$RBANS_attent_FU12 <- NULL
nap_2LP_pass$RBANS_total_FU12 <- NULL
nap_2LP_pass$total_TAU_FU12 <- NULL
nap_2LP_pass$phospho_TAU_FU12 <- NULL
nap_2LP_pass$b_amyloid_FU12 <- NULL
nap_2LP_pass_long = nap_2LP_pass[,c(1, 7, 10, 15, 17, 20, 23, 31, 32, 33, 34, 35, 36, 42, 44, 46, 48, 50, 52, 54, 56, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69)]

#### GM CBF
## organize data (grep GM CBF only and convert table to a longitudinal one)
NAP_2LP_pass_GM_CBF = nap_2LP_pass_long[,c(1,2,3,4,5,6,7,14,18)]
GM_CBF_long = melt(NAP_2LP_pass_GM_CBF, id=c(1,2,3,4,5,6,7), value.name="GM_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_GM <- ggplot(data=GM_CBF_long, aes(x= TimePoint, y= GM_Avg_CBF, group=DCCID))
p_GM + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM Precuneus 
NAP_2LP_pass_Prec_CBF = nap_2LP_pass_long[,c(1,2,3,4,5,6,7,15,19)]
Prec_CBF_long = melt(NAP_2LP_pass_Prec_CBF, id=c(1,2,3,4,5,6,7), value.name="Prec_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_Prec <- ggplot(data=Prec_CBF_long, aes(x= TimePoint, y=Prec_Avg_CBF, group=DCCID))
p_Prec + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM PCC 
NAP_2LP_pass_PCC_CBF = nap_2LP_pass_long[,c(1,2,3,4,5,6,7,16,20)]
PCC_CBF_long = melt(NAP_2LP_pass_PCC_CBF, id=c(1,2,3,4,5,6,7), value.name="PCC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_PCC <- ggplot(data=PCC_CBF_long, aes(x= TimePoint, y=PCC_Avg_CBF, group=DCCID))
p_PCC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)

#### GM ACC
NAP_2LP_pass_ACC_CBF = nap_2LP_pass_long[,c(1,2,3,4,5,6,7,17,21)]
ACC_CBF_long = melt(NAP_2LP_pass_ACC_CBF, id=c(1,2,3,4,5,6,7), value.name="ACC_Avg_CBF", variable.name="TimePoint")
## set the base for the spaghetti plot
p_ACC <- ggplot(data=ACC_CBF_long, aes(x= TimePoint, y=ACC_Avg_CBF, group=DCCID))
p_ACC + geom_line() + stat_smooth(aes(group=1)) + stat_summary(aes(group=1), geom="point", fun.y=mean, shape = 17, size =3) + facet_grid(. ~ E4_allele_no)



# create table with up till 24 month cohort pass CBF data only
cohort_FU24_pass = subset(data, Project=="PRE" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU12=="pass" & Precuneus_QC_FU12=="pass" & PCC_QC_FU12=="pass" & ACC_QC_FU12=="pass" & GM_QC_FU24=="pass" & Precuneus_QC_FU24=="pass" & PCC_QC_FU24=="pass" & ACC_QC_FU24=="pass")

# create table with up till 12 month cohort pass CBF data only
cohort_FU12_pass = subset(data, Project=="PRE" & GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass" & GM_QC_FU12=="pass" & Precuneus_QC_FU12=="pass" & PCC_QC_FU12=="pass" & ACC_QC_FU12=="pass")

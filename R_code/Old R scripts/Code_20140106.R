# clear workspace
rm(list = ls())

# load modules
require(foreign)
require(MASS)

# load datasets
data	= read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/QC/for analyses/ALL_BL_data_20150119.csv",h=T,sep=",")
volumes = read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/QC/volumes spreadsheets/PreventAD_volumes_BL_20141127.csv", h=T, sep=",")
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




# create tables with pass BL CBF data only
# GM_BL_pass	= subset(data, GM_QC_BL=="pass")
# Prec_BL_pass= subset(data, Precuneus_QC_BL=="pass")
# PCC_BL_pass	= subset(data, PCC_QC_BL=="pass")
# ACC_BL_pass	= subset(data, ACC_QC_BL=="pass")
ASL_BL_pass = subset(data, GM_QC_BL=="pass" & Precuneus_QC_BL=="pass" & PCC_QC_BL=="pass" & ACC_QC_BL=="pass")
# Add column with age in years 
ASL_BL_pass$Age_BL_years = ASL_BL_pass$Age_BL_scan_date/12

### Cleanup data to keep only BL stuff
ASL_BL_pass$ACC_Avg_FU03 <- NULL
ASL_BL_pass$ACC_QC_FU03 <- NULL
ASL_BL_pass$PCC_Avg_FU03 <- NULL
ASL_BL_pass$PCC_QC_FU03 <- NULL
ASL_BL_pass$Precuneus_Avg_FU03 <- NULL
ASL_BL_pass$Precuneus_QC_FU03 <- NULL
ASL_BL_pass$GM_Avg_FU03 <- NULL
ASL_BL_pass$GM_QC_FU03 <- NULL
ASL_BL_pass$ACC_Avg_FU24 <- NULL
ASL_BL_pass$ACC_QC_FU24 <- NULL
ASL_BL_pass$PCC_Avg_FU24 <- NULL
ASL_BL_pass$PCC_QC_FU24 <- NULL
ASL_BL_pass$Precuneus_Avg_FU24 <- NULL
ASL_BL_pass$Precuneus_QC_FU24 <- NULL
ASL_BL_pass$GM_Avg_FU24 <- NULL
ASL_BL_pass$GM_QC_FU24 <- NULL
ASL_BL_pass$ACC_Avg_FU12 <- NULL
ASL_BL_pass$ACC_QC_FU12 <- NULL
ASL_BL_pass$PCC_Avg_FU12 <- NULL
ASL_BL_pass$PCC_QC_FU12 <- NULL
ASL_BL_pass$Precuneus_Avg_FU12 <- NULL
ASL_BL_pass$Precuneus_QC_FU12 <- NULL
ASL_BL_pass$GM_Avg_FU12 <- NULL
ASL_BL_pass$GM_QC_FU12 <- NULL
ASL_BL_pass$RBANS_immed_mem_FU03 <- NULL
ASL_BL_pass$RBANS_visuo_const_FU03 <- NULL
ASL_BL_pass$RBANS_language_FU03 <- NULL
ASL_BL_pass$RBANS_attent_FU03 <- NULL
ASL_BL_pass$RBANS_delay_mem_FU03 <- NULL
ASL_BL_pass$RBANS_total_FU03 <- NULL
ASL_BL_pass$RBANS_immed_mem_FU12 <- NULL
ASL_BL_pass$RBANS_visuo_const_FU12 <- NULL
ASL_BL_pass$RBANS_language_FU12 <- NULL
ASL_BL_pass$RBANS_attent_FU12 <- NULL
ASL_BL_pass$RBANS_delay_mem_FU12 <- NULL
ASL_BL_pass$RBANS_total_FU12 <- NULL
ASL_BL_pass$total_TAU_FU03 <- NULL
ASL_BL_pass$phospho_TAU_FU03 <- NULL
ASL_BL_pass$b_amyloid_FU03 <- NULL
ASL_BL_pass$total_TAU_FU12 <- NULL
ASL_BL_pass$phospho_TAU_FU12 <- NULL
ASL_BL_pass$b_amyloid_FU12 <- NULL


##### AGE EFFECT #####
## 1. look at relation between age and Mean GM CBF
par(mfrow=c(2,2))
model = ASL_BL_pass$GM_Avg_BL~ASL_BL_pass$Age_BL_years
plot(model, ylim=c(30,80), xlim=c(55,80), xlab="Age (years)", ylab="Grey Matter Average CBF (ml/100g/min)", col="red", font.lab=4, cex.lab=1.3)
fit = rlm(model)
abline(fit, col="blue")
summary(fit)
# compute pvalue
# p_value = 2*min(1-pnorm(Value/Std), pnorm(Value/Std)) where Value and Std can be found in summary(fit)
# Value/Std = t value of the regression table
# p_value = 2*min(1-pnorm(tval), pnorm(tval))
model2 = ASL_BL_pass$Precuneus_Avg_BL~ASL_BL_pass$Age_BL_years
plot(model2, ylim=c(30,80), xlim=c(55,80), xlab="Age (years)", ylab="Precuneus Average CBF (ml/100g/min)", col="red", font.lab=4, cex.lab=1.3)
fit2 = rlm(model2)
abline(fit2, col="blue")
summary(fit2)
model3 = ASL_BL_pass$PCC_Avg_BL~ASL_BL_pass$Age_BL_years
plot(model3, ylim=c(30,80), xlim=c(55,80), xlab="Age (years)", ylab="PCC Average CBF (ml/100g/min)", col="red", font.lab=4, cex.lab=1.3)
fit3 = rlm(model3)
abline(fit3, col="blue")
summary(fit3)
model4 = ASL_BL_pass$ACC_Avg_BL~ASL_BL_pass$Age_BL_years
plot(model4, ylim=c(30,80), xlim=c(55,80), xlab="Age (years)", ylab="ACC Average CBF (ml/100g/min)", col="red", font.lab=4, cex.lab=1.3)
fit4 = rlm(model4)
abline(fit4, col="blue")
summary(fit4)

## 2. relation between age, gender and mean CBF
female 	= subset(ASL_BL_pass, Gender=="Female")
male 	= subset(ASL_BL_pass, Gender=="Male")

# A. For the graph
par(mfrow=c(2,2))
# GM
plot(female$GM_Avg_BL~female$Age_BL_years, col="red", ylim=c(20,100))
points(male$GM_Avg_BL~male$Age_BL_years, col="blue")
model.female	= female$GM_Avg_BL~female$Age_BL_years
model.male		= male$GM_Avg_BL~male$Age_BL_years
fit.female		= rlm(model.female)
fit.male		= rlm(model.male)
abline(fit.female, col="red")
abline(fit.male, col="blue")
legend(55,30, c("female", "male"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# Precuneus
plot(female$Precuneus_Avg_BL~female$Age_BL_years, col="red", ylim=c(20,100))
points(male$Precuneus_Avg_BL~male$Age_BL_years, col="blue")
model.female	= female$Precuneus_Avg_BL~female$Age_BL_years
model.male		= male$Precuneus_Avg_BL~male$Age_BL_years
fit.female		= rlm(model.female)
fit.male		= rlm(model.male)
abline(fit.female, col="red")
abline(fit.male, col="blue")
legend(55,30, c("female", "male"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# PCC
plot(female$PCC_Avg_BL~female$Age_BL_years, col="red", ylim=c(20,100))
points(male$PCC_Avg_BL~male$Age_BL_years, col="blue")
model.female	= female$PCC_Avg_BL~female$Age_BL_years
model.male		= male$PCC_Avg_BL~male$Age_BL_years
fit.female		= rlm(model.female)
fit.male		= rlm(model.male)
abline(fit.female, col="red")
abline(fit.male, col="blue")
legend(55,30, c("female", "male"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# ACC
plot(female$ACC_Avg_BL~female$Age_BL_years, col="red", ylim=c(20,100))
points(male$ACC_Avg_BL~male$Age_BL_years, col="blue")
model.female	= female$ACC_Avg_BL~female$Age_BL_years
model.male		= male$ACC_Avg_BL~male$Age_BL_years
fit.female		= rlm(model.female)
fit.male		= rlm(model.male)
abline(fit.female, col="red")
abline(fit.male, col="blue")
legend(55,30, c("female", "male"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))

# B. models
# GM
model_GM = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$Gender_bin
fit_GM	 = rlm(model_GM)
summary(fit_GM)
pvalue_gm_age    = 2*min(1-pnorm(-3.2252), pnorm(-3.2252))
pvalue_gm_gender = 2*min(1-pnorm(0.5283), pnorm(0.5283))
pvalue_gm_inter  = 2*min(1-pnorm(-0.7203), pnorm(-0.7203))
# Precuneus
model_Prec = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$Gender_bin
fit_Prec   = rlm(model_Prec)
summary(fit_Prec)
pvalue_prec_age    = 2*min(1-pnorm(-3.5024), pnorm(-3.5024))
pvalue_prec_gender = 2*min(1-pnorm(-0.9688), pnorm(-0.9688))
pvalue_prec_inter  = 2*min(1-pnorm(0.7473), pnorm(0.7473))
# PCC
model_PCC = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$Gender_bin
fit_PCC   = rlm(model_PCC)
summary(fit_PCC)
pvalue_pcc_age    = 2*min(1-pnorm(-3.6693), pnorm(-3.6693))
pvalue_pcc_gender = 2*min(1-pnorm(0.3389), pnorm(0.3389))
pvalue_pcc_inter  = 2*min(1-pnorm(-0.6633), pnorm(-0.6633))
# PCC
model_ACC = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$Gender_bin
fit_ACC   = rlm(model_ACC)
summary(fit_ACC)
pvalue_acc_age    = 2*min(1-pnorm(-3.0425), pnorm(-3.0425))
pvalue_acc_gender = 2*min(1-pnorm(0.2765), pnorm(0.2765))
pvalue_acc_inter  = 2*min(1-pnorm(-0.2555), pnorm(-0.2555))

## 3. Relation between genotype and CBF
E4pos   = subset(ASL_BL_pass, E4_allele_no==1)
E4neg   = subset(ASL_BL_pass, E4_allele_no==0)
BDNFpos = subset(ASL_BL_pass, BDNF_copie_no==1)
BDNFneg = subset(ASL_BL_pass, BDNF_copie_no==0)

# A. ApoE graphs
par(mfrow=c(2,2))
# GM
plot(E4pos$GM_Avg_BL~E4pos$Age_BL_years, col="red", ylim=c(20,100))
points(E4neg$GM_Avg_BL~ E4neg$Age_BL_years, col="blue")
model.E4pos	= E4pos$GM_Avg_BL~E4pos$Age_BL_years
model.E4neg	= E4neg$GM_Avg_BL~E4neg$Age_BL_years
fit.E4pos	= rlm(model.E4pos)
fit.E4neg	= rlm(model.E4neg)
abline(fit.E4pos, col="red")
abline(fit.E4neg, col="blue")
legend(56,100, c("ApoE+", "ApoE-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# Precuneus
plot(E4pos$Precuneus_Avg_BL~E4pos$Age_BL_years, col="red", ylim=c(20,100))
points(E4neg$Precuneus_Avg_BL~E4neg$Age_BL_years, col="blue")
model.E4pos	= E4pos$Precuneus_Avg_BL~E4pos$Age_BL_years
model.E4neg		= E4neg$Precuneus_Avg_BL~E4neg$Age_BL_years
fit.E4pos		= rlm(model.E4pos)
fit.E4neg		= rlm(model.E4neg)
abline(fit.E4pos, col="red")
abline(fit.E4neg, col="blue")
legend(56,100, c("ApoE+", "ApoE-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# PCC
plot(E4pos$PCC_Avg_BL~E4pos$Age_BL_years, col="red", ylim=c(20,100))
points(E4neg$PCC_Avg_BL~E4neg$Age_BL_years, col="blue")
model.E4pos	= E4pos$PCC_Avg_BL~E4pos$Age_BL_years
model.E4neg		= E4neg$PCC_Avg_BL~E4neg$Age_BL_years
fit.E4pos		= rlm(model.E4pos)
fit.E4neg		= rlm(model.E4neg)
abline(fit.E4pos, col="red")
abline(fit.male, col="blue")
legend(56,100, c("ApoE+", "ApoE-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# ACC
plot(E4pos$ACC_Avg_BL~E4pos$Age_BL_years, col="red", ylim=c(20,100))
points(E4neg$ACC_Avg_BL~E4neg$Age_BL_years, col="blue")
model.E4pos	= E4pos$ACC_Avg_BL~E4pos$Age_BL_years
model.E4neg	= E4neg$ACC_Avg_BL~E4neg$Age_BL_years
fit.E4pos	= rlm(model.E4pos)
fit.E4neg	= rlm(model.E4neg)
abline(fit.E4pos, col="red")
abline(fit.E4neg, col="blue")
legend(56,100, c("ApoE+", "ApoE-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))

# B. ApoE models
# GM
model_GM_apo = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$E4_allele_no
fit_GM_apo	 = rlm(model_GM_apo)
summary(fit_GM_apo)
pvalue_gm_age    = 2*min(1-pnorm(-4.0279), pnorm(-4.0279))
pvalue_gm_apo = 2*min(1-pnorm(-1.1512), pnorm(-1.1512))
pvalue_gm_inter  = 2*min(1-pnorm(1.2704), pnorm(1.2704))
# Precuneus
model_Prec_apo = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$E4_allele_no
fit_Prec_apo   = rlm(model_Prec_apo)
summary(fit_Prec_apo)
pvalue_prec_age    = 2*min(1-pnorm(-3.5605), pnorm(-3.5605))
pvalue_prec_apo = 2*min(1-pnorm(0.6407), pnorm(0.6407))
pvalue_prec_inter  = 2*min(1-pnorm(-0.7191), pnorm(-0.7191))
# PCC
model_PCC_apo = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$E4_allele_no
fit_PCC_apo   = rlm(model_PCC_apo)
summary(fit_PCC_apo)
pvalue_pcc_age    = 2*min(1-pnorm(-4.2465), pnorm(-4.2465))
pvalue_pcc_apo = 2*min(1-pnorm(-1.4281), pnorm(-1.4281))
pvalue_pcc_inter  = 2*min(1-pnorm(1.4871), pnorm(1.4871))
# PCC
model_ACC_apo = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$E4_allele_no
fit_ACC_apo   = rlm(model_ACC_apo)
summary(fit_ACC_apo)
pvalue_acc_age    = 2*min(1-pnorm(-3.3560), pnorm(-3.3560))
pvalue_acc_apo = 2*min(1-pnorm(-0.8504), pnorm(-0.8504))
pvalue_acc_inter  = 2*min(1-pnorm(0.8905), pnorm(0.8905))

# Removing the candidates with protective Intron M factor (Intron M negative is protective)
ASL_BL_pass_IntronMpos = subset(ASL_BL_pass, Intron_M_Protective_Variant ==0)
# GM
model_GM_apo = ASL_BL_pass_IntronMpos$GM_Avg_BL ~ ASL_BL_pass_IntronMpos$Age_BL_years * ASL_BL_pass_IntronMpos$E4_allele_no
fit_GM_apo	 = rlm(model_GM_apo)
summary(fit_GM_apo)
pvalue_gm_age    = 2*min(1-pnorm(-2.6149), pnorm(-2.6149))
pvalue_gm_apo = 2*min(1-pnorm(-0.9440), pnorm(-0.9440))
pvalue_gm_inter  = 2*min(1-pnorm(1.0357), pnorm(1.0357))
# Precuneus
model_Prec_apo = ASL_BL_pass_IntronMpos$Precuneus_Avg_BL ~ ASL_BL_pass_IntronMpos$Age_BL_years * ASL_BL_pass_IntronMpos$E4_allele_no
fit_Prec_apo   = rlm(model_Prec_apo)
summary(fit_Prec_apo)
pvalue_prec_age    = 2*min(1-pnorm(-2.2548), pnorm(-2.2548))
pvalue_prec_apo = 2*min(1-pnorm(-0.5285), pnorm(-0.5285))
pvalue_prec_inter  = 2*min(1-pnorm(0.5484), pnorm(0.5484))
# PCC
model_PCC_apo = ASL_BL_pass_IntronMpos$PCC_Avg_BL ~ ASL_BL_pass_IntronMpos$Age_BL_years * ASL_BL_pass_IntronMpos$E4_allele_no
fit_PCC_apo   = rlm(model_PCC_apo)
summary(fit_PCC_apo)
pvalue_pcc_age    = 2*min(1-pnorm(-2.4468), pnorm(-2.4468))
pvalue_pcc_apo = 2*min(1-pnorm(-0.8907), pnorm(-0.8907))
pvalue_pcc_inter  = 2*min(1-pnorm(0.9197), pnorm(0.9197))
# PCC
model_ACC_apo = ASL_BL_pass_IntronMpos$ACC_Avg_BL ~ ASL_BL_pass_IntronMpos$Age_BL_years * ASL_BL_pass_IntronMpos$E4_allele_no
fit_ACC_apo   = rlm(model_ACC_apo)
summary(fit_ACC_apo)
pvalue_acc_age    = 2*min(1-pnorm(-2.9396), pnorm(-2.9396))
pvalue_acc_apo = 2*min(1-pnorm(-0.9539), pnorm(-0.9539))
pvalue_acc_inter  = 2*min(1-pnorm(0.9789), pnorm(0.9789))

# Looking at differences between E4+ K+ and E4+ K-
E4pos_Kpos= subset(E4pos, K_variant_positive==1)
E4pos_Kneg= subset(E4pos, K_variant_positive==0)
wilcox.test(E4pos_Kpos$GM_Avg_BL,E4pos_Kneg$GM_Avg_BL)
wilcox.test(E4pos_Kpos$Precuneus_Avg_BL,E4pos_Kneg$Precuneus_Avg_BL)
wilcox.test(E4pos_Kpos$PCC_Avg_BL,E4pos_Kneg$PCC_Avg_BL)
wilcox.test(E4pos_Kpos$ACC_Avg_BL,E4pos_Kneg$ACC_Avg_BL)
# graphs
# boxplots
boxplot(E4pos$GM_Avg_BL~E4pos$K_variant_positive, main="Grey Matter", ylab="Avg CBF", xlab="K variant positive")
boxplot(E4pos$Precuneus_Avg_BL~E4pos$K_variant_positive, main="Precuneus", ylab="Avg CBF", xlab="K variant positive")
boxplot(E4pos$PCC_Avg_BL~E4pos$K_variant_positive, main="PCC", ylab="Avg CBF", xlab="K variant positive")
boxplot(E4pos$ACC_Avg_BL~E4pos$K_variant_positive, main="ACC", ylab="Avg CBF", xlab="K variant positive")
# scatterplots
plot(E4pos_Kpos$GM_Avg_BL~E4pos_Kpos$Age_BL_years,col="red", ylim=c(20,100))
points(E4pos_Kneg$GM_Avg_BL~E4pos_Kneg$Age_BL_years,col="blue")
legend(56,100, c("ApoE+K+", "ApoE+K-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
model.EposKpos_GM = rlm(E4pos_Kpos$GM_Avg_BL~E4pos_Kpos$Age_BL_years)
model.EposKneg_GM = rlm(E4pos_Kneg$GM_Avg_BL~E4pos_Kneg$Age_BL_years)
abline(model.EposKpos_GM, col="red")
abline(model.EposKneg_GM, col="blue")
plot(E4pos_Kpos$Precuneus_Avg_BL~E4pos_Kpos$Age_BL_years,col="red", ylim=c(20,100))
points(E4pos_Kneg$Precuneus_Avg_BL~E4pos_Kneg$Age_BL_years,col="blue")
model.EposKpos_Prec = rlm(E4pos_Kpos$Precuneus_Avg_BL~E4pos_Kpos$Age_BL_years)
model.EposKneg_Prec = rlm(E4pos_Kneg$Precuneus_Avg_BL~E4pos_Kneg$Age_BL_years)
abline(model.EposKpos_Prec, col="red")
abline(model.EposKneg_Prec, col="blue")
legend(56,100, c("ApoE+K+", "ApoE+K-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
plot(E4pos_Kpos$PCC_Avg_BL~E4pos_Kpos$Age_BL_years,col="red", ylim=c(20,100))
points(E4pos_Kneg$PCC_Avg_BL~E4pos_Kneg$Age_BL_years,col="blue")
model.EposKpos_PCC = rlm(E4pos_Kpos$PCC_Avg_BL~E4pos_Kpos$Age_BL_years)
model.EposKneg_PCC = rlm(E4pos_Kneg$PCC_Avg_BL~E4pos_Kneg$Age_BL_years)
abline(model.EposKpos_PCC, col="red")
abline(model.EposKneg_PCC, col="blue")
legend(56,100, c("ApoE+K+", "ApoE+K-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
plot(E4pos_Kpos$ACC_Avg_BL~E4pos_Kpos$Age_BL_years,col="red", ylim=c(20,100))
points(E4pos_Kneg$ACC_Avg_BL~E4pos_Kneg$Age_BL_years,col="blue")
model.EposKpos_ACC = rlm(E4pos_Kpos$PCC_Avg_BL~E4pos_Kpos$Age_BL_years)
model.EposKneg_ACC = rlm(E4pos_Kneg$PCC_Avg_BL~E4pos_Kneg$Age_BL_years)
abline(model.EposKpos_ACC, col="red")
abline(model.EposKneg_ACC, col="blue")
legend(56,100, c("ApoE+K+", "ApoE+K-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# models
# GM
model_GM_E4K = E4pos$GM_Avg_BL ~ E4pos$Age_BL_years * E4pos$K_variant_positive 
fit_GM_E4K = rlm(model_GM_E4K)
summary(fit_GM_E4K)
pvalue_gm_intercept = 2*min(1-pnorm(2.7243), pnorm(2.7243))
pvalue_gm_age       = 2*min(1-pnorm(-0.0436), pnorm(-0.0436))
pvalue_gm_e4k 	    = 2*min(1-pnorm(1.0336), pnorm(1.0336))
pvalue_gm_inter     = 2*min(1-pnorm(-0.9535), pnorm(-0.9535))
# Precuneus
model_Prec_E4K = E4pos$Precuneus_Avg_BL ~ E4pos$Age_BL_years * E4pos$K_variant_positive 
fit_Prec_E4K = rlm(model_Prec_E4K)
summary(fit_Prec_E4K)
pvalue_prec_intercept = 2*min(1-pnorm(2.8932), pnorm(2.8932))
pvalue_prec_age    	  = 2*min(1-pnorm(-0.8039), pnorm(-0.8039))
pvalue_prec_e4k 	  = 2*min(1-pnorm(0.4330), pnorm(0.4330))
pvalue_prec_inter 	  = 2*min(1-pnorm(-0.3407), pnorm(-0.3407))
# PCC
model_PCC_E4K = E4pos$PCC_Avg_BL ~ E4pos$Age_BL_years * E4pos$K_variant_positive 
fit_PCC_E4K = rlm(model_PCC_E4K)
summary(fit_PCC_E4K)
pvalue_pcc_intercept = 2*min(1-pnorm(1.6585), pnorm(1.6585))
pvalue_pcc_age   	 = 2*min(1-pnorm(0.2682), pnorm(0.2682))
pvalue_pcc_e4k 	  	 = 2*min(1-pnorm(1.5321), pnorm(1.5321))
pvalue_pcc_inter 	 = 2*min(1-pnorm(-1.4476), pnorm(-1.4476))
# ACC
model_ACC_E4K = E4pos$ACC_Avg_BL ~ E4pos$Age_BL_years * E4pos$K_variant_positive 
fit_ACC_E4K = rlm(model_ACC_E4K)
summary(fit_ACC_E4K)
pvalue_pcc_intercept = 2*min(1-pnorm(2.6609), pnorm(2.6609))
pvalue_acc_age    	 = 2*min(1-pnorm(0.1301), pnorm(0.1301))
pvalue_acc_e4k 	 	 = 2*min(1-pnorm(0.9326), pnorm(0.9326))
pvalue_acc_inter 	 = 2*min(1-pnorm(-0.8443), pnorm(-0.8443))



# C. BDNF graphs
par(mfrow=c(2,2))
# GM
plot(BDNFpos$GM_Avg_BL~BDNFpos$Age_BL_years, col="red", ylim=c(20,100))
points(BDNFneg$GM_Avg_BL~ BDNFneg$Age_BL_years, col="blue")
model.BDNFpos	= BDNFpos$GM_Avg_BL~BDNFpos$Age_BL_years
model.BDNFneg	= BDNFneg$GM_Avg_BL~BDNFneg$Age_BL_years
fit.BDNFpos	= rlm(model.BDNFpos)
fit.BDNFneg	= rlm(model.BDNFneg)
abline(fit.BDNFpos, col="red")
abline(fit.BDNFneg, col="blue")
legend(70,100, c("BDNF+", "BDNF-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# Precuneus
plot(BDNFpos$Precuneus_Avg_BL~BDNFpos$Age_BL_years, col="red", ylim=c(20,100))
points(BDNFneg$Precuneus_Avg_BL~BDNFneg$Age_BL_years, col="blue")
model.BDNFpos	= BDNFpos$Precuneus_Avg_BL~BDNFpos$Age_BL_years
model.BDNFneg		= BDNFneg$Precuneus_Avg_BL~BDNFneg$Age_BL_years
fit.BDNFpos		= rlm(model.BDNFpos)
fit.BDNFneg		= rlm(model.BDNFneg)
abline(fit.BDNFpos, col="red")
abline(fit.BDNFneg, col="blue")
legend(70,100, c("BDNF+", "BDNF-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# PCC
plot(BDNFpos$PCC_Avg_BL~BDNFpos$Age_BL_years, col="red", ylim=c(20,100))
points(BDNFneg$PCC_Avg_BL~BDNFneg$Age_BL_years, col="blue")
model.BDNFpos	= BDNFpos$PCC_Avg_BL~BDNFpos$Age_BL_years
model.BDNFneg		= BDNFneg$PCC_Avg_BL~BDNFneg$Age_BL_years
fit.BDNFpos		= rlm(model.BDNFpos)
fit.BDNFneg		= rlm(model.BDNFneg)
abline(fit.BDNFpos, col="red")
abline(fit.male, col="blue")
legend(70,100, c("BDNF+", "BDNF-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
# ACC
plot(BDNFpos$ACC_Avg_BL~BDNFpos$Age_BL_years, col="red", ylim=c(20,100))
points(BDNFneg$ACC_Avg_BL~BDNFneg$Age_BL_years, col="blue")
model.BDNFpos	= BDNFpos$ACC_Avg_BL~BDNFpos$Age_BL_years
model.BDNFneg	= BDNFneg$ACC_Avg_BL~BDNFneg$Age_BL_years
fit.BDNFpos	= rlm(model.BDNFpos)
fit.BDNFneg	= rlm(model.BDNFneg)
abline(fit.BDNFpos, col="red")
abline(fit.BDNFneg, col="blue")
legend(70,100, c("BDNF+", "BDNF-"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))

# D. BDNF models
model_GM_bdnf = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$BDNF_copie_no
fit_GM_bdnf	 = rlm(model_GM_bdnf)
summary(fit_GM_bdnf)
pvalue_gm_age    = 2*min(1-pnorm(-2.7299), pnorm(-2.7299))
pvalue_gm_bdnf   = 2*min(1-pnorm(0.8726), pnorm(0.8726))
pvalue_gm_inter  = 2*min(1-pnorm(-0.9419), pnorm(-0.9419))
# Precuneus
model_Prec_bdnf = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$BDNF_copie_no
fit_Prec_bdnf   = rlm(model_Prec_bdnf)
summary(fit_Prec_bdnf)
pvalue_prec_age    = 2*min(1-pnorm(-2.2048), pnorm(-2.2048))
pvalue_prec_bdnf   = 2*min(1-pnorm(1.5554), pnorm(1.5554))
pvalue_prec_inter  = 2*min(1-pnorm(-1.6243), pnorm(-1.6243))
# PCC
model_PCC_bdnf = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$BDNF_copie_no
fit_PCC_bdnf   = rlm(model_PCC_bdnf)
summary(fit_PCC_bdnf)
pvalue_pcc_age    = 2*min(1-pnorm(-2.5818), pnorm(-2.5818))
pvalue_pcc_bdnf   = 2*min(1-pnorm(0.9855), pnorm(0.9855))
pvalue_pcc_inter  = 2*min(1-pnorm(-1.0849), pnorm(-1.0849))
# PCC
model_ACC_bdnf = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$Age_BL_years * ASL_BL_pass$BDNF_copie_no
fit_ACC_bdnf   = rlm(model_ACC_bdnf)
summary(fit_ACC_bdnf)
pvalue_acc_age    = 2*min(1-pnorm(-2.3608), pnorm(-2.3608))
pvalue_acc_bdnf   = 2*min(1-pnorm(0.5367), pnorm(0.5367))
pvalue_acc_inter  = 2*min(1-pnorm(-0.6522), pnorm(-0.6522))


# E. ApoI_BL
# divide variable in high and low ApoI
median(ASL_BL_pass$ApoAI_BL, na.rm=T)
ApoAIpos = subset(ASL_BL_pass, ApoAI_BL > median(ASL_BL_pass$ApoAI_BL, na.rm=T))
ApoAIneg = subset(ASL_BL_pass, ApoAI_BL < median(ASL_BL_pass$ApoAI_BL, na.rm=T))

# GM model
model_GM_ApoAI  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$ApoAI_BL *ASL_BL_pass$Age_BL_years
fit_GM_ApoAI 	= rlm(model_GM_ApoAI)
summary(fit_GM_ApoAI)
# Precuneus model
model_Prec_ApoAI  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$ApoAI_BL *ASL_BL_pass$Age_BL_years
fit_Prec_ApoAI 	= rlm(model_Prec_ApoAI)
summary(fit_Prec_ApoAI)
# PCC model
model_PCC_ApoAI  = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$ApoAI_BL *ASL_BL_pass$Age_BL_years
fit_PCC_ApoAI 	= rlm(model_PCC_ApoAI)
summary(fit_PCC_ApoAI)
#ACC model
model_ACC_ApoAI  = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$ApoAI_BL *ASL_BL_pass$Age_BL_years
fit_ACC_ApoAI 	= rlm(model_ACC_ApoAI)
summary(fit_ACC_ApoAI)
# GM graph
plot(ApoAIpos$GM_Avg_BL~ApoAIpos$Age_BL_years, col="red")
points(ApoAIneg$GM_Avg_BL~ApoAIneg$Age_BL_years, col="blue")
model_ApoAIpos = ApoAIpos$GM_Avg_BL~ApoAIpos$Age_BL_years
model_ApoAIneg = ApoAIneg$GM_Avg_BL~ApoAIneg$Age_BL_years
fit_ApoAIpos = rlm(model_ApoAIpos)
fit_ApoAIneg = rlm(model_ApoAIneg)
abline(fit_ApoAIpos, col="red")
abline(fit_ApoAIneg, col="blue")


# F. RBANS
# models
# total score 
model_GM_total  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_total_BL *ASL_BL_pass$Age_BL_years
fit_GM_total 	= rlm(model_GM_total)
summary(fit_GM_total)
model_Prec_total  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_total_BL *ASL_BL_pass$Age_BL_years
fit_Prec_total 	= rlm(model_Prec_total)
summary(fit_Prec_total)
model_PCC_total  = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$RBANS_total_BL *ASL_BL_pass$Age_BL_years
fit_PCC_total 	= rlm(model_PCC_total)
summary(fit_PCC_total)
model_ACC_total  = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$RBANS_total_BL *ASL_BL_pass$Age_BL_years
fit_ACC_total 	= rlm(model_ACC_total)
summary(fit_ACC_total)
# immed_memory
model_GM_immed_mem  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_immed_mem_BL *ASL_BL_pass$Age_BL_years
fit_GM_immed_mem 	= rlm(model_GM_immed_mem)
summary(fit_GM_immed_mem)
model_Prec_immed_mem  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_immed_mem_BL *ASL_BL_pass$Age_BL_years
fit_Prec_immed_mem 	= rlm(model_Prec_immed_mem)
summary(fit_Prec_immed_mem)
model_PCC_immed_mem  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_immed_mem_BL *ASL_BL_pass$Age_BL_years
fit_PCC_immed_mem 	= rlm(model_PCC_immed_mem)
summary(fit_PCC_immed_mem)
model_ACC_immed_mem  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_immed_mem_BL *ASL_BL_pass$Age_BL_years
fit_ACC_immed_mem 	= rlm(model_ACC_immed_mem)
summary(fit_ACC_immed_mem)
# attention
model_GM_attent  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_attent_BL *ASL_BL_pass$Age_BL_years
fit_GM_attent 	= rlm(model_GM_attent)
summary(fit_GM_attent)
model_Prec_attent  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_attent_BL *ASL_BL_pass$Age_BL_years
fit_Prec_attent 	= rlm(model_Prec_attent)
summary(fit_Prec_attent)
model_PCC_attent  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_attent_BL *ASL_BL_pass$Age_BL_years
fit_PCC_attent 	= rlm(model_PCC_attent)
summary(fit_PCC_attent)
model_ACC_attent  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_attent_BL *ASL_BL_pass$Age_BL_years
fit_ACC_attent 	= rlm(model_ACC_attent)
summary(fit_ACC_attent)
# language
model_GM_lang  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_language_BL *ASL_BL_pass$Age_BL_years
fit_GM_lang 	= rlm(model_GM_lang)
summary(fit_GM_lang)
model_Prec_lang  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_language_BL *ASL_BL_pass$Age_BL_years
fit_Prec_lang 	= rlm(model_Prec_lang)
summary(fit_Prec_lang)
model_PCC_lang  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_language_BL *ASL_BL_pass$Age_BL_years
fit_PCC_lang 	= rlm(model_PCC_lang)
summary(fit_PCC_lang)
model_ACC_lang  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_language_BL *ASL_BL_pass$Age_BL_years
fit_ACC_lang 	= rlm(model_ACC_lang)
summary(fit_ACC_lang)
# visuo_const
model_GM_visuo_const  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_visuo_const_BL *ASL_BL_pass$Age_BL_years
fit_GM_visuo_const 	= rlm(model_GM_visuo_const)
summary(fit_GM_visuo_const)
model_Prec_visuo_const  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_visuo_const_BL *ASL_BL_pass$Age_BL_years
fit_Prec_visuo_const 	= rlm(model_Prec_visuo_const)
summary(fit_Prec_visuo_const)
model_PCC_visuo_const  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_visuo_const_BL *ASL_BL_pass$Age_BL_years
fit_PCC_visuo_const 	= rlm(model_PCC_visuo_const)
summary(fit_PCC_visuo_const)
model_ACC_visuo_const  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_visuo_const_BL *ASL_BL_pass$Age_BL_years
fit_ACC_visuo_const 	= rlm(model_ACC_visuo_const)
summary(fit_ACC_visuo_const)

# delayed mem
model_GM_delay_mem  = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$RBANS_delay_mem_BL *ASL_BL_pass$Age_BL_years
fit_GM_delay_mem 	= rlm(model_GM_delay_mem)
summary(fit_GM_delay_mem)
pvalue_del_mem_BL = 2*min(1-pnorm(2.3494), pnorm(2.3494))
pvalue_age 		  = 2*min(1-pnorm(1.8502), pnorm(1.8502))
pvalue_inter	  = 2*min(1-pnorm(-2.1323), pnorm(-2.1323))

model_Prec_delay_mem  = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$RBANS_delay_mem_BL *ASL_BL_pass$Age_BL_years
fit_Prec_delay_mem 	= rlm(model_Prec_delay_mem)
summary(fit_Prec_delay_mem)
pvalue_del_mem_BL = 2*min(1-pnorm(3.1455), pnorm(3.1455))
pvalue_age	      = 2*min(1-pnorm(2.7306), pnorm(2.7306))
pvalue_inter      = 2*min(1-pnorm(-2.9625), pnorm(-2.9625))

model_PCC_delay_mem  = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$RBANS_delay_mem_BL *ASL_BL_pass$Age_BL_years
fit_PCC_delay_mem 	= rlm(model_PCC_delay_mem)
summary(fit_PCC_delay_mem)
pvalue_del_mem_BL = 2*min(1-pnorm(2.1449), pnorm(2.1449))
pvalue_age	      = 2*min(1-pnorm(1.7348), pnorm(1.7348))
pvalue_inter      = 2*min(1-pnorm(-1.9706), pnorm(-1.9706))

model_ACC_delay_mem  = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$RBANS_delay_mem_BL *ASL_BL_pass$Age_BL_years
fit_ACC_delay_mem 	= rlm(model_ACC_delay_mem)
summary(fit_ACC_delay_mem)
pvalue_del_mem_BL = 2*min(1-pnorm(2.6067), pnorm(2.6067))
pvalue_age	      = 2*min(1-pnorm(2.3517), pnorm(2.3517))
pvalue_inter      = 2*min(1-pnorm(-2.5189), pnorm(-2.5189))


### Graphs for delayed memory
par(mfrow=c(2,2))
high_delay_score = subset(ASL_BL_pass, RBANS_delay_mem_BL > median(ASL_BL_pass$RBANS_delay_mem_BL, na.rm=T))
low_delay_score = subset(ASL_BL_pass, RBANS_delay_mem_BL < median(ASL_BL_pass$RBANS_delay_mem_BL, na.rm=T))

plot(high_delay_score$GM_Avg_BL~high_delay_score$Age_BL_years, col="red", ylim=c(20,110))
points(low_delay_score $GM_Avg_BL~ low_delay_score $Age_BL_years, col="blue")
model_high_delay = high_delay_score$GM_Avg_BL~high_delay_score$Age_BL_years
model_low_delay = low_delay_score $GM_Avg_BL~low_delay_score $Age_BL_years
fit_high_delay = rlm(model_high_delay)
fit_low_delay = rlm(model_low_delay)
abline(fit_high_delay, col="red")
abline(fit_low_delay, col="blue")
legend(63,110, c("High DM Score", "Low DM Score"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
plot(high_delay_score$Precuneus_Avg_BL~high_delay_score$Age_BL_years, col="red", ylim=c(20,110))
points(low_delay_score $Precuneus_Avg_BL~ low_delay_score $Age_BL_years, col="blue")
model_high_delay = high_delay_score$Precuneus_Avg_BL~high_delay_score$Age_BL_years
model_low_delay = low_delay_score $Precuneus_Avg_BL~low_delay_score $Age_BL_years
fit_high_delay = rlm(model_high_delay)
fit_low_delay = rlm(model_low_delay)
abline(fit_high_delay, col="red")
abline(fit_low_delay, col="blue")
legend(63,110, c("High DM Score", "Low DM Score"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
plot(high_delay_score$PCC_Avg_BL~high_delay_score$Age_BL_years, col="red", ylim=c(20,110))
points(low_delay_score $PCC_Avg_BL~ low_delay_score $Age_BL_years, col="blue")
model_high_delay = high_delay_score$PCC_Avg_BL~high_delay_score$Age_BL_years
model_low_delay = low_delay_score $PCC_Avg_BL~low_delay_score $Age_BL_years
fit_high_delay = rlm(model_high_delay)
fit_low_delay = rlm(model_low_delay)
abline(fit_high_delay, col="red")
abline(fit_low_delay, col="blue")
legend(63,110, c("High DM Score", "Low DM Score"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))
plot(high_delay_score$ACC_Avg_BL~high_delay_score$Age_BL_years, col="red", ylim=c(20,110))
points(low_delay_score $ACC_Avg_BL~ low_delay_score $Age_BL_years, col="blue")
model_high_delay = high_delay_score$ACC_Avg_BL~high_delay_score$Age_BL_years
model_low_delay = low_delay_score $ACC_Avg_BL~low_delay_score $Age_BL_years
fit_high_delay = rlm(model_high_delay)
fit_low_delay = rlm(model_low_delay)
abline(fit_high_delay, col="red")
abline(fit_low_delay, col="blue")
legend(63,110, c("High DM Score", "Low DM Score"), lty=c(1,1), lwd=c(1.5,1.5), col=c("red","blue"))

### G. HC volumes
volumes$HC_left_scaled     = volumes$hippocampus_left * volumes$ScaleFactor
volumes$HC_right_scaled    = volumes$hippocampus_right * volumes$ScaleFactor
volumes$HC_combined_scaled = volumes$HC_left_scaled + volumes$HC_right_scaled
# Append volumes to ASL_BL_pass
HC			 = volumes[,c(1,46,47,48)]
colnames(HC) = c("DCCID","HC_right_scaled", "HC_left_scaled", "HC_combined_scaled")
ASL_BL_pass  = merge(ASL_BL_pass, HC, by="DCCID")

# Modeling
model_GM_HC = ASL_BL_pass$GM_Avg_BL ~ ASL_BL_pass$HC_combined_scaled * ASL_BL_pass$Age_BL_years
fit_GM_HC = rlm(model_GM_HC)
summary(fit_GM_HC)
model_Prec_HC = ASL_BL_pass$Precuneus_Avg_BL ~ ASL_BL_pass$HC_combined_scaled * ASL_BL_pass$Age_BL_years
fit_Prec_HC = rlm(model_Prec_HC)
summary(fit_Prec_HC)
model_PCC_HC = ASL_BL_pass$PCC_Avg_BL ~ ASL_BL_pass$HC_combined_scaled * ASL_BL_pass$Age_BL_years
fit_PCC_HC = rlm(model_PCC_HC)
summary(fit_PCC_HC)
model_ACC_HC = ASL_BL_pass$ACC_Avg_BL ~ ASL_BL_pass$HC_combined_scaled * ASL_BL_pass$Age_BL_years
fit_ACC_HC = rlm(model_ACC_HC)
summary(fit_ACC_HC)
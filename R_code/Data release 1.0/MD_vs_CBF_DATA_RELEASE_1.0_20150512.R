###### SOURCE LoadDataRelease.R FILE and DATA_RELEASE_1.0_20150512
source('~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/R_code/LoadDataRelease.R', chdir = TRUE)

###### LOAD SPREADSHEETS
LoadSpreadsheets( '~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/for analyses/DATA_RELEASE_1.0', 
				  '~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/QC spreadsheets/DKT_template/QC_spreadsheet_20150508.csv',
				  '~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/for analyses/DKT template/new_MD_CBF_DKT_results/ASL_MD_OHBM_dataset.csv'
				)

# load medication data (the spreadsheet below has been appended the changes noticed at baseline)
LoadAdditionalSpreadsheet('curr_med_use_data','~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/for analyses/DATA_RELEASE_1.0_spreadsheets_BLonly_passPreCBF/current_med_use_with_changes_at_baseline.csv')




# shows available dataframe that has been cleaned-up
ls( pattern = "data" )

#### LOOK AT DATA THAT ARE PASSED FOR BOTH MD AND CBF DATA

# Check normality of CBF and MD data
hist( MD_data$asl_avg_GM )
hist( MD_data$asl_PCN )
hist( MD_data$MD_pcn)
	# => everyting has a normal distribution
	
## General_physical
cbf_general_physical = merge( general_physical_data, MD_data, by="CandID" )
hist( cbf_general_physical$bmi )
# check relation between BMI and average CBF in GM
plot(asl_PCN ~ bmi,cbf_general_physical)
abline(lm(asl_PCN ~ bmi,cbf_general_physical))
summary(lm(asl_PCN ~ bmi,cbf_general_physical))
# check relation between BMI and average CBF in PCN
plot(asl_PCN ~ bmi,cbf_general_physical)
abline(lm(asl_PCN ~ bmi,cbf_general_physical))
summary(lm(asl_PCN ~ bmi,cbf_general_physical))
	# => no relation between CBF and BMI

## General_medical_history
cbf_general_med_hist = merge( general_medical_history_data, MD_data, by="CandID" )
# hypertension treatment vs CBF
plot(asl_avg_GM ~ X3_treatment_hypertension, cbf_general_med_hist)
plot(asl_PCN ~ X3_treatment_hypertension, cbf_general_med_hist)
	# => no effect on CBF 
# hyperlipidemia vs CBF
plot(asl_avg_GM ~ X4_treatment_hyperlipidemia, cbf_general_med_hist)
plot(asl_PCN ~ X4_treatment_hyperlipidemia, cbf_general_med_hist)
	# => no effect on CBF 
# diabetes
plot(asl_avg_GM ~ X5_treatment_diabetes, cbf_general_med_hist)
plot(asl_PCN ~ X5_treatment_diabetes, cbf_general_med_hist)
	# => no effect on CBF
# tia
plot(asl_avg_GM ~ X7_tia, cbf_general_med_hist)
	# => only 4 subjects with tia, no effect on CBF
# head injury
plot(asl_avg_GM ~ X8_head_injury, cbf_general_med_hist)
	# => no effect on CBF
# depression
cbf_general_med_hist$X10_b_depression_specify_bin = ifelse( is.na(cbf_general_med_hist$X10_b_depression_specify), "no", "yes")
boxplot(asl_avg_GM ~X10_b_depression_specify_bin , cbf_general_med_hist)
	# => no effect on CBF
# myocardial infraction
cbf_general_med_hist$X10_e_myocardial_infraction_specify_bin = ifelse( is.na(cbf_general_med_hist$X10_e_myocardial_infraction_specify), "no", "yes")
boxplot(asl_avg_GM ~ X10_e_myocardial_infraction_specify_bin , cbf_general_med_hist)
	# => no effect on CBF

## general_information_data
cbf_general_info = merge( general_information_data, MD_data, by="CandID" )
# sleeping habits
plot(asl_avg_GM ~ X12_sleeping_habits, cbf_general_info)
	# => no effect on CBF
# cigarette smoking habits
plot(asl_avg_GM ~ X13_cigarette_smoking, cbf_general_info)
plot(asl_PCN ~ X13_cigarette_smoking, cbf_general_info)
	# => no effect on CBF
# number of years of education 	
plot(asl_avg_GM ~ X4_education, cbf_general_info)
abline(lm(asl_avg_GM ~ X4_education, cbf_general_info))
plot(asl_PCN ~ X4_education, cbf_general_info)
abline(lm(asl_PCN ~ X4_education, cbf_general_info))
	=> no effect on CBF
# level of education
plot(asl_avg_GM ~ level_of_education, cbf_general_info)
plot(asl_PCN ~ level_of_education, cbf_general_info)
	=> no effect on CBF
# exercice levels
plot(asl_avg_GM ~ X22_activity_level, cbf_general_info)
plot(asl_PCN ~ X22_activity_level, cbf_general_info)
	=> maybe not enough subjects to evaluate this?
	
## cohort_candidate_info_data
cbf_cohort_candidate_info_data = merge( cohort_candidate_info_data, MD_data, by="CandID" )
# gender
plot(asl_avg_GM ~ Gender, cbf_cohort_candidate_info_data)
cbf_cohort_candidate_info_data_males   = subset(cbf_cohort_candidate_info_data, Gender == "Male")
cbf_cohort_candidate_info_data_females = subset(cbf_cohort_candidate_info_data, Gender == "Female")
hist(cbf_cohort_candidate_info_data_males$asl_avg_GM)
hist(cbf_cohort_candidate_info_data_females$asl_avg_GM)
hist(cbf_cohort_candidate_info_data_males$asl_PCN)
hist(cbf_cohort_candidate_info_data_females$asl_PCN)
	# => males distribution is sort of normal...
t.test(asl_avg_GM ~ Gender, cbf_cohort_candidate_info_data)
wilcox.test(asl_avg_GM ~ Gender, cbf_cohort_candidate_info_data)
	# => significant effect of Gender on CBF
plot(asl_avg_GM ~ age_EN, cbf_cohort_candidate_info_data_males , col="blue", xlab="Candidate's age (years)", ylab="Grey matter average CBF (ml/100g/min)")
points(asl_avg_GM ~ age_EN, cbf_cohort_candidate_info_data_females , col="red")
abline(lm(asl_avg_GM ~ age_EN, cbf_cohort_candidate_info_data_males), col="blue")
abline(lm(asl_avg_GM ~ age_EN, cbf_cohort_candidate_info_data_females), col="red")
legend(55, 40, c("females", "males"), col=c("red", "blue"), pch="o")
summary(lm(asl_avg_GM ~ age_EN * Gender, cbf_cohort_candidate_info_data)) 
summary(lm(asl_avg_GM ~ age_EN + Gender, cbf_cohort_candidate_info_data)) 
plot(asl_PCN ~ age_EN, cbf_cohort_candidate_info_data_males , col="blue", xlab="Candidate's age (years)", ylab="Average CBF in Precuneus (ml/100g/min)")
points(asl_PCN ~ age_EN, cbf_cohort_candidate_info_data_females , col="red")
abline(lm(asl_PCN ~ age_EN, cbf_cohort_candidate_info_data_males), col="blue")
abline(lm(asl_PCN ~ age_EN, cbf_cohort_candidate_info_data_females), col="red")
legend(55, 40, c("females", "males"), col=c("red", "blue"), pch="o")
summary(lm(asl_PCN ~ age_EN * Gender, cbf_cohort_candidate_info_data)) 
summary(lm(asl_PCN ~ age_EN + Gender, cbf_cohort_candidate_info_data)) 
plot(MD_pcn ~ age_EN, cbf_cohort_candidate_info_data_males , col="blue", xlab="Candidate's age (years)", ylab="MD in Precuneus")
points(MD_pcn ~ age_EN, cbf_cohort_candidate_info_data_females , col="red")
abline(lm(MD_pcn ~ age_EN, cbf_cohort_candidate_info_data_males), col="blue")
abline(lm(MD_pcn ~ age_EN, cbf_cohort_candidate_info_data_females), col="red")
legend(55, 40, c("females", "males"), col=c("red", "blue"), pch="o")
summary(lm(MD_pcn ~ age_EN * Gender, cbf_cohort_candidate_info_data)) 
summary(lm(MD_pcn ~ age_EN + Gender, cbf_cohort_candidate_info_data)) 

## family_history_data
cbf_family_history_data = merge( family_history_data, MD_data, by="CandID" )
cbf_family_history_data$age_from_family_onset = cbf_family_history_data$family_mean_age_onset - cbf_family_history_data$age_EN
plot(asl_avg_GM ~ age_from_family_onset, cbf_family_history_data)
boxplot(asl_avg_GM ~ number_of_family_member_with_AD, cbf_family_history_data)
	# => no effect of age from family onset on CBF or number of members in the family with AD on CBF
	
## AD8_data
AD8_data_BL  = subset(AD8_data, Visit_label %in% c("NAPBL00","PREBL00"))
cbf_AD8_data = merge( AD8_data_BL, MD_data, by="CandID") 
boxplot(asl_avg_GM~total_score, cbf_AD8_data)
	# => no effect on CBF
	
## baseline_data (blood pressure data)
cbf_baseline_data = merge( baseline_data, MD_data, by="CandID")
plot(asl_avg_GM ~ bp_systolic, cbf_baseline_data, xlab="Blood pressure systolic (mm Hg)", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ bp_systolic, cbf_baseline_data))
summary(lm(asl_avg_GM ~ bp_systolic, cbf_baseline_data))
plot(asl_avg_GM ~ bp_diastolic, cbf_baseline_data, xlab="Blood pressure diastolic (mm Hg)", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ bp_diastolic, cbf_baseline_data))
summary(lm(asl_avg_GM ~ bp_diastolic, cbf_baseline_data))
plot(asl_avg_GM ~ X9_pulse, cbf_baseline_data, xlab="Pulse", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X9_pulse, cbf_baseline_data))
summary(lm(asl_avg_GM ~ X9_pulse, cbf_baseline_data))
plot(asl_avg_GM ~ X8_bp_normal, cbf_baseline_data)
plot(asl_avg_GM ~ X9_pulse_regularity, cbf_baseline_data)
	# => negative correlation between bp systolic and GM CBF (trend)
	# => negative correlation between bp diastolic and GM CBF (trend)
	# => positive correlation between pulse and GM CBF
plot(asl_PCN ~ bp_systolic, cbf_baseline_data, xlab="Blood pressure systolic (mm Hg)", ylab="Average CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ bp_systolic, cbf_baseline_data))
summary(lm(asl_PCN ~ bp_systolic, cbf_baseline_data))
plot(asl_PCN ~ bp_diastolic, cbf_baseline_data, xlab="Blood pressure diastolic (mm Hg)", ylab="Average CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ bp_diastolic, cbf_baseline_data))
summary(lm(asl_PCN ~ bp_diastolic, cbf_baseline_data))
plot(asl_PCN ~ X9_pulse, cbf_baseline_data, xlab="Pulse", ylab="Average CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X9_pulse, cbf_baseline_data))
summary(lm(asl_PCN ~ X9_pulse, cbf_baseline_data))

# enrollment_data (EKG + blood pressure)
cbf_enrollment_data = merge( enrollment_data, MD_data, by="CandID")
# replace not_answered by NA
cbf_enrollment_data$X6_EKG_result[cbf_enrollment_data$X6_EKG_result == "not_answered"] = NA
factor(cbf_enrollment_data$X6_EKG_result)
plot(asl_avg_GM ~ X6_EKG_result, cbf_enrollment_data, xlab="EKG result", ylab="Grey matter average CBF (ml/100g/min)")
t.test(cbf_enrollment_data$asl_avg_GM ~ cbf_enrollment_data$X6_EKG_result)		# => significant
wilcox.test(cbf_enrollment_data$asl_avg_GM ~ cbf_enrollment_data$X6_EKG_result)  # => trend
cbf_enrollment_data$bp_systolic  = as.character(cbf_enrollment_data$bp_systolic )
cbf_enrollment_data$bp_diastolic = as.character(cbf_enrollment_data$bp_diastolic)
cbf_enrollment_data$bp_systolic  = as.numeric(cbf_enrollment_data$bp_systolic )
cbf_enrollment_data$bp_diastolic = as.numeric(cbf_enrollment_data$bp_diastolic)
plot(asl_avg_GM ~ bp_systolic, cbf_enrollment_data, xlab="Blood pressure systolic (mm Hg)", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ bp_systolic, cbf_enrollment_data))
summary(lm(asl_avg_GM ~ bp_systolic, cbf_enrollment_data))
plot(asl_avg_GM ~ bp_diastolic, cbf_enrollment_data, xlab="Blood pressure diastolic (mm Hg)", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ bp_diastolic, cbf_enrollment_data))
summary(lm(asl_avg_GM ~ bp_diastolic, cbf_enrollment_data))
plot(asl_avg_GM ~ X9_pulse, cbf_enrollment_data, xlab="Pulse", ylab="Grey matter average CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X9_pulse, cbf_enrollment_data))
summary(lm(asl_avg_GM ~ X9_pulse, cbf_enrollment_data))
# Remove outlier and plot again graphs above
cbf_enrollment_data[144,] = NA
	# => negative correlation between bp systolic and GM CBF 
	# => negative correlation between bp diastolic and GM CBF (trend)
	# => no correlation between pulse and GM CBF

## genetics_data
cbf_genetics_data = merge(genetics_data, MD_data, by="CandID")
# create binary variables
#ApoE_bin = ifelse( grepl( "4" , cbf_genetics_data$ApoE 			  	  ), 1,  0 )
#BchE_bin = ifelse( grepl( "A" , cbf_genetics_data$BcheE_K_variant 	  ), 1,  0 )
#BDNF_bin = ifelse( grepl( "A" , cbf_genetics_data$BDNF 			  	  ), 1,  0 )
#TLR4_bin = ifelse( grepl( "G" , cbf_genetics_data$TLR4_rs_4986790 	  ), 1,  0 )
#PPP2_bin = ifelse( grepl( "TT", cbf_genetics_data$PPP2r1A_rs_10406151 ), 1,  0 )
# check that already existing variable are correct
#summary( ApoE_bin == cbf_genetics_data$E4_allele_Bin    )
#summary( BchE_bin == cbf_genetics_data$K_variant_bin    )
#summary( BDNF_bin == cbf_genetics_data$BDNF_copie_bin   )
#summary( TLR4_bin == cbf_genetics_data$TLR4_allele_no   )
#summary( PPP2_bin == cbf_genetics_data$ppp2r1A_copie_no )
plot       ( asl_avg_GM ~ E4_allele_Bin, cbf_genetics_data )
t.test     ( asl_avg_GM ~ E4_allele_Bin, cbf_genetics_data )
wilcox.test( asl_avg_GM ~ E4_allele_Bin, cbf_genetics_data )
plot       ( asl_avg_GM ~ ApoE_HMGR_bin, cbf_genetics_data )
t.test     ( asl_avg_GM ~ ApoE_HMGR_bin, cbf_genetics_data )
wilcox.test( asl_avg_GM ~ ApoE_HMGR_bin, cbf_genetics_data )

## handedness_data
cbf_handedness_data = merge(handedness_data, MD_data, by="CandID")
plot(asl_avg_GM ~ interpretation, cbf_handedness_data)

## lab_data
lab_data_BL  = subset(lab_data,    Visit_label %in% c("NAPBL00","PREBL00"))
cbf_lab_data = merge( lab_data_BL, MD_data, by="CandID") 
plot(asl_avg_GM ~ X1_hb_result, cbf_lab_data, xlab="Hemoglobin result (g/L)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X1_hb_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X1_hb_result, cbf_lab_data))
plot(asl_PCN ~ X1_hb_result, cbf_lab_data, xlab="Hemoglobin result (g/L)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X1_hb_result, cbf_lab_data))
summary(lm(asl_PCN ~ X1_hb_result, cbf_lab_data))
	# => strong effect of HB levels on GM CBF & PCN CBF
plot(asl_avg_GM ~ X2_ht_result, cbf_lab_data, xlab="Hematocrit", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X2_ht_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X2_ht_result, cbf_lab_data))
plot(asl_PCN ~ X2_ht_result, cbf_lab_data, xlab="Hematocrit", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X2_ht_result, cbf_lab_data))
summary(lm(asl_PCN ~ X2_ht_result, cbf_lab_data))
	# => strong effect of HT levels on GM CBF & PCN CBF
plot(asl_avg_GM ~ X3_wbc_result, cbf_lab_data, xlab="White blood cell (x10_9)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X3_wbc_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X3_wbc_result, cbf_lab_data))
plot(asl_PCN ~ X3_wbc_result, cbf_lab_data, xlab="White blood cell (x10_9)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X3_wbc_result, cbf_lab_data))
summary(lm(asl_PCN ~ X3_wbc_result, cbf_lab_data))
	# => trend of WBC on GM CBF but no effet on PCN CBF
plot(asl_avg_GM ~ X4_rbc_result, cbf_lab_data, xlab="White blood cell (x10_12)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X4_rbc_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X4_rbc_result, cbf_lab_data))
plot(asl_PCN ~ X4_rbc_result, cbf_lab_data, xlab="White blood cell (x10_12)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X4_rbc_result, cbf_lab_data))
summary(lm(asl_PCN ~ X4_rbc_result, cbf_lab_data))
	# => strong effect of RBC on GM CBF & PCN CBF
plot(asl_avg_GM ~ X5_platelet_result, cbf_lab_data, xlab="Platelets (x10_9)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X5_platelet_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X5_platelet_result, cbf_lab_data))
plot(asl_PCN ~ X5_platelet_result, cbf_lab_data, xlab="Platelets (x10_9)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X5_platelet_result, cbf_lab_data))
summary(lm(asl_PCN ~ X5_platelet_result, cbf_lab_data))
	# => small effect of platelets on GM CBF and PCN CBF
plot(asl_avg_GM ~ X6_sodium_result, cbf_lab_data, xlab="Sodium (millimole)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X6_sodium_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X6_sodium_result, cbf_lab_data))
plot(asl_PCN ~ X6_sodium_result, cbf_lab_data, xlab="Sodium (millimole)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X6_sodium_result, cbf_lab_data))
summary(lm(asl_PCN ~ X6_sodium_result, cbf_lab_data))
	# => no effect of sodium on CBF
cbf_lab_data$X7_potassium_result = as.numeric(cbf_lab_data$X7_potassium_result)
plot(asl_avg_GM ~ X7_potassium_result, cbf_lab_data, xlab="Potassium (millimole)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X7_potassium_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X7_potassium_result, cbf_lab_data))
plot(asl_PCN ~ X7_potassium_result, cbf_lab_data, xlab="Potassium (millimole)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X7_potassium_result, cbf_lab_data))
summary(lm(asl_PCN ~ X7_potassium_result, cbf_lab_data))
	# => no effect of potassium on CBF
plot(asl_avg_GM ~ X8_chloride_result, cbf_lab_data, xlab="Chloride (millimole)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X8_chloride_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X8_chloride_result, cbf_lab_data))
plot(asl_PCN ~ X8_chloride_result, cbf_lab_data, xlab="Chloride (millimole)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X8_chloride_result, cbf_lab_data))
summary(lm(asl_PCN ~ X8_chloride_result, cbf_lab_data))
	# => no effect of chloride on CBF
plot(asl_avg_GM ~ X9_urea_result, cbf_lab_data, xlab="Urea (millimole)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X9_urea_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X9_urea_result, cbf_lab_data))
plot(asl_PCN ~ X9_urea_result, cbf_lab_data, xlab="Urea (millimole)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X9_urea_result, cbf_lab_data))
summary(lm(asl_PCN ~ X9_urea_result, cbf_lab_data))
	# => no effect of urea on CBF
plot(asl_avg_GM ~ X10_creatinine_result, cbf_lab_data, xlab="Creatinine (millimole)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X10_creatinine_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X10_creatinine_result, cbf_lab_data))
plot(asl_PCN ~ X10_creatinine_result, cbf_lab_data, xlab="Creatinine (millimole)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X10_creatinine_result, cbf_lab_data))
summary(lm(asl_PCN ~ X10_creatinine_result, cbf_lab_data))
	# => trend between creatinine and CBF (GM & PCN)
plot(asl_avg_GM ~ X11_protein_tot_result, cbf_lab_data, xlab="Total protein (grams)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X11_protein_tot_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X11_protein_tot_result, cbf_lab_data))
plot(asl_PCN ~ X11_protein_tot_result, cbf_lab_data, xlab="Total protein (grams)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X11_protein_tot_result, cbf_lab_data))
summary(lm(asl_PCN ~ X11_protein_tot_result, cbf_lab_data))
	# => no effect of total protein on CBF
cbf_lab_data$X12_bilirubin_tot_result = as.numeric(cbf_lab_data$X12_bilirubin_tot_result)
plot(asl_avg_GM ~ X12_bilirubin_tot_result, cbf_lab_data, xlab="Total bilirubin (micromoles)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X12_bilirubin_tot_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X12_bilirubin_tot_result, cbf_lab_data))
plot(asl_PCN ~ X12_bilirubin_tot_result, cbf_lab_data, xlab="Total bilirubin (micromoles)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X12_bilirubin_tot_result, cbf_lab_data))
summary(lm(asl_PCN ~ X12_bilirubin_tot_result, cbf_lab_data))
	# => no effect of bilirubin on CBF
plot(asl_avg_GM ~ X13_alk_phosph_result, cbf_lab_data, xlab="Alkaline phosphatase", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X13_alk_phosph_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X13_alk_phosph_result, cbf_lab_data))
plot(asl_PCN ~ X13_alk_phosph_result, cbf_lab_data, xlab="Alkaline phosphatase", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X13_alk_phosph_result, cbf_lab_data))
summary(lm(asl_PCN ~ X13_alk_phosph_result, cbf_lab_data))
	# => no effect of alkaline phosphatase on CBF
plot(asl_avg_GM ~ X14_alt_result, cbf_lab_data, xlab="Alanine aminotransferase", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X14_alt_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X14_alt_result, cbf_lab_data))
plot(asl_PCN ~ X14_alt_result, cbf_lab_data, xlab="Alanine aminotransferase", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X14_alt_result, cbf_lab_data))
summary(lm(asl_PCN ~ X14_alt_result, cbf_lab_data))
	# => trend between alt and PCN CBF (but nothing with GM_CBF)
plot(asl_avg_GM ~ X15_ast_result, cbf_lab_data, xlab="Aspartate aminotransferase", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X15_ast_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X15_ast_result, cbf_lab_data))
plot(asl_PCN ~ X15_ast_result, cbf_lab_data, xlab="Aspartate aminotransferase", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X15_ast_result, cbf_lab_data))
summary(lm(asl_PCN ~ X15_ast_result, cbf_lab_data))
	# => no effect of ast on CBF
plot(asl_avg_GM ~ X16_albumin_result, cbf_lab_data, xlab="Albumin (grams)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X16_albumin_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X16_albumin_result, cbf_lab_data))
plot(asl_PCN ~ X16_albumin_result, cbf_lab_data, xlab="Albumin (grams)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X16_albumin_result, cbf_lab_data))
summary(lm(asl_PCN ~ X16_albumin_result, cbf_lab_data))
	=> trend of albumin on CBF
plot(asl_avg_GM ~ X17_ggt_result, cbf_lab_data, xlab="Gamma-glutamyl transferase", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X17_ggt_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X17_ggt_result, cbf_lab_data))
plot(asl_PCN ~ X17_ggt_result, cbf_lab_data, xlab="Gamma-glutamyl transferase", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X17_ggt_result, cbf_lab_data))
summary(lm(asl_PCN ~ X17_ggt_result, cbf_lab_data))
	# => no effect of GGT on CBF
cbf_lab_data$X18_hba1c_result = as.numeric(cbf_lab_data$X18_hba1c_result)
plot(asl_avg_GM ~ X18_hba1c_result, cbf_lab_data, xlab="Hemoglobin A1c (grams)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X18_hba1c_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X18_hba1c_result, cbf_lab_data))
plot(asl_PCN ~ X18_hba1c_result, cbf_lab_data, xlab="Hemoglobin A1c (grams)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X18_hba1c_result, cbf_lab_data))
summary(lm(asl_PCN ~ X18_hba1c_result, cbf_lab_data))
	# => small effect of hba1c on CBF
plot(asl_avg_GM ~ X19_tsh_result, cbf_lab_data, xlab="Thyroid-stimulating hormone (milliunits)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X19_tsh_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X19_tsh_result, cbf_lab_data))
plot(asl_PCN ~ X19_tsh_result, cbf_lab_data, xlab="Thyroid-stimulating hormone (milliunits)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X19_tsh_result, cbf_lab_data))
summary(lm(asl_PCN ~ X19_tsh_result, cbf_lab_data))
 	# => effect of TSH on CBF
plot(asl_avg_GM ~ X20_pt_result, cbf_lab_data, xlab="Prothrombin time (seconds)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X20_pt_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X20_pt_result, cbf_lab_data))
plot(asl_PCN ~ X20_pt_result, cbf_lab_data, xlab="Prothrombin time (seconds)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X20_pt_result, cbf_lab_data))
summary(lm(asl_PCN ~ X20_pt_result, cbf_lab_data))
 	# => no effect of PT on CBF
plot(asl_avg_GM ~ X21_inr_result, cbf_lab_data, xlab="International normalized ratio", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X21_inr_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X21_inr_result, cbf_lab_data))
plot(asl_PCN ~ X21_inr_result, cbf_lab_data, xlab="International normalized ratio", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X21_inr_result, cbf_lab_data))
summary(lm(asl_PCN ~ X21_inr_result, cbf_lab_data))
	# => no effect of INR on CBF
plot(asl_avg_GM ~ X22_ptt_result, cbf_lab_data, xlab="Partial thromboplastin time (seconds)", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ X22_ptt_result, cbf_lab_data))
summary(lm(asl_avg_GM ~ X22_ptt_result, cbf_lab_data))
plot(asl_PCN ~ X22_ptt_result, cbf_lab_data, xlab="Partial thromboplastin time (seconds)", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ X22_ptt_result, cbf_lab_data))
summary(lm(asl_PCN ~ X22_ptt_result, cbf_lab_data))
	# => no effect of PTT on CBF

## LP_data
LP_data_BL  = subset(LP_data,    Visit_label %in% c("NAPLP00","PRELP00"))
cbf_LP_data = merge( LP_data_BL, MD_data, by="CandID") 
# Create BP systolic & diastolic variables from bp_sitting
cbf_LP_data$bp_sitting = as.character( cbf_LP_data$bp_sitting)
cbf_LP_data$bp_systolic_sitting  = as.numeric(sapply(strsplit(cbf_LP_data$bp_sitting, "/"), function(x) x[1]))
cbf_LP_data$bp_diastolic_sitting = as.numeric(sapply(strsplit(cbf_LP_data$bp_sitting, "/"), function(x) x[2]))
plot(asl_avg_GM ~ bp_systolic_sitting, cbf_LP_data)
abline(lm(asl_avg_GM ~ bp_systolic_sitting, cbf_LP_data))
summary(lm(asl_avg_GM ~ bp_systolic_sitting, cbf_LP_data))
plot(asl_PCN ~ bp_systolic_sitting, cbf_LP_data)
abline(lm(asl_PCN ~ bp_systolic_sitting, cbf_LP_data))
summary(lm(asl_PCN ~ bp_systolic_sitting, cbf_LP_data))

plot(asl_avg_GM ~ pulse_sitting, cbf_LP_data)
abline(lm(asl_avg_GM ~ pulse_sitting, cbf_LP_data))
	# => no effect of pulse on CBF
plot(asl_avg_GM ~ X2_rbc_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X2_rbc_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X2_rbc_result, cbf_LP_data))
plot(asl_PCN ~ X2_rbc_result, cbf_LP_data)
abline(lm(asl_PCN ~ X2_rbc_result, cbf_LP_data))
summary(lm(asl_PCN ~ X2_rbc_result, cbf_LP_data))
	# => no effect of RBC on CBF
plot(asl_avg_GM ~ X3_neutrophils_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X3_neutrophils_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X3_neutrophils_result, cbf_LP_data))
plot(asl_PCN ~ X3_neutrophils_result, cbf_LP_data)
abline(lm(asl_PCN ~ X3_neutrophils_result, cbf_LP_data))
summary(lm(asl_PCN ~ X3_neutrophils_result, cbf_LP_data))
	# => no effect of neutrophils on CBF
plot(asl_avg_GM ~ X4_lymphocytes_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X4_lymphocytes_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X4_lymphocytes_result, cbf_LP_data))
plot(asl_PCN ~ X4_lymphocytes_result, cbf_LP_data)
abline(lm(asl_PCN ~ X4_lymphocytes_result, cbf_LP_data))
summary(lm(asl_PCN ~ X4_lymphocytes_result, cbf_LP_data))
	# => no effect of lymphocytes on CBF
plot(asl_avg_GM ~ X6_epithelial_cells_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6_epithelial_cells_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6_epithelial_cells_result, cbf_LP_data))
plot(asl_PCN ~ X6_epithelial_cells_result, cbf_LP_data)
abline(lm(asl_PCN ~ X6_epithelial_cells_result, cbf_LP_data))
summary(lm(asl_PCN ~ X6_epithelial_cells_result, cbf_LP_data))
	# => no effect on epithelial cells on CBF
plot(asl_avg_GM ~ X8_wbc_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X8_wbc_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X8_wbc_result, cbf_LP_data))
plot(asl_PCN ~ X8_wbc_result, cbf_LP_data)
abline(lm(asl_PCN ~ X8_wbc_result, cbf_LP_data))
summary(lm(asl_PCN ~ X8_wbc_result, cbf_LP_data))
	# => no effect on WBC on CBF
plot(asl_avg_GM ~ ELISA_tau, cbf_LP_data)
abline(lm(asl_avg_GM ~ ELISA_tau, cbf_LP_data))
summary(lm(asl_avg_GM ~ ELISA_tau, cbf_LP_data))
plot(asl_PCN ~ ELISA_tau, cbf_LP_data)
abline(lm(asl_PCN ~ ELISA_tau, cbf_LP_data))
summary(lm(asl_PCN ~ ELISA_tau, cbf_LP_data))
	# => no effect of tau on CBF
plot(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data)
abline(lm(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data))
summary(lm(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data))
plot(asl_PCN ~ ELISA_b_amyloid, cbf_LP_data)
abline(lm(asl_PCN ~ ELISA_b_amyloid, cbf_LP_data))
summary(lm(asl_PCN ~ ELISA_b_amyloid, cbf_LP_data))
	# => no effect of b_amyloid on CBF
plot(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data)
abline(lm(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data))
summary(lm(asl_avg_GM ~ ELISA_b_amyloid, cbf_LP_data))
plot(asl_PCN ~ ELISA_ptau, cbf_LP_data)
abline(lm(asl_PCN ~ ELISA_ptau, cbf_LP_data))
summary(lm(asl_PCN ~ ELISA_ptau, cbf_LP_data))
	# => no effect of ptau on CBF
plot(asl_avg_GM ~ tau_over_abeta_ratio, cbf_LP_data)
abline(lm(asl_avg_GM ~ tau_over_abeta_ratio, cbf_LP_data))
summary(lm(asl_avg_GM ~ tau_over_abeta_ratio, cbf_LP_data))
plot(asl_PCN ~ tau_over_abeta_ratio, cbf_LP_data)
abline(lm(asl_PCN ~ tau_over_abeta_ratio, cbf_LP_data))
summary(lm(asl_PCN ~ tau_over_abeta_ratio, cbf_LP_data))
	# no effect of tau/abeta on CBF
plot(asl_avg_GM ~ ELISA_NFL, cbf_LP_data)
abline(lm(asl_avg_GM ~ ELISA_NFL, cbf_LP_data))
summary(lm(asl_avg_GM ~ ELISA_NFL, cbf_LP_data))
plot(asl_PCN ~ ELISA_NFL, cbf_LP_data)
abline(lm(asl_PCN ~ ELISA_NFL, cbf_LP_data))
summary(lm(asl_PCN ~ ELISA_NFL, cbf_LP_data))
	# => no effect of NFL on CBF
plot(asl_avg_GM ~ X6plex_ApoAI, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoAI, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoAI, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoAI, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoAI, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoAI, cbf_LP_data))
	# => no effect of ApoAI on CBF
plot(asl_avg_GM ~ X6plex_ApoAII, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoAII, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoAII, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoAII, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoAII, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoAII, cbf_LP_data))
	# => no effect of ApoAII on CBF
plot(asl_avg_GM ~ X6plex_ApoB, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoB, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoB, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoB, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoB, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoB, cbf_LP_data))
	# => no effect of ApoB on CBF
plot(asl_avg_GM ~ X6plex_ApoCII, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoCII, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoCII, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoCII, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoCII, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoCII, cbf_LP_data))
	# => no effect of ApoCII on CBF
plot(asl_avg_GM ~ X6plex_ApoCIII, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoCIII, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoCIII, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoCIII, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoCIII, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoCIII, cbf_LP_data))
	# => no effect of ApoCIII on CBF
plot(asl_avg_GM ~ X6plex_ApoE, cbf_LP_data)
abline(lm(asl_avg_GM ~ X6plex_ApoE, cbf_LP_data))
summary(lm(asl_avg_GM ~ X6plex_ApoE, cbf_LP_data))
plot(asl_PCN ~ X6plex_ApoE, cbf_LP_data)
abline(lm(asl_PCN ~ X6plex_ApoE, cbf_LP_data))
summary(lm(asl_PCN ~ X6plex_ApoE, cbf_LP_data))
	# => no effect of ApoE on CBF
	## ===> probably need more subjects though
plot(asl_avg_GM ~ X9_glucose_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X9_glucose_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X9_glucose_result, cbf_LP_data))
plot(asl_PCN ~ X9_glucose_result, cbf_LP_data)
abline(lm(asl_PCN ~ X9_glucose_result, cbf_LP_data))
summary(lm(asl_PCN ~ X9_glucose_result, cbf_LP_data))
	# => no effect of glucose on CBF
plot(asl_avg_GM ~ X10_micro_protein_result, cbf_LP_data)
abline(lm(asl_avg_GM ~ X10_micro_protein_result, cbf_LP_data))
summary(lm(asl_avg_GM ~ X10_micro_protein_result, cbf_LP_data))
plot(asl_PCN ~ X10_micro_protein_result, cbf_LP_data)
abline(lm(asl_PCN ~ X10_micro_protein_result, cbf_LP_data))
summary(lm(asl_PCN ~ X10_micro_protein_result, cbf_LP_data))

## RBANS_data
RBANS_data_BL  = subset(RBANS_data,    Visit_label %in% c("NAPBL00","PREBL00"))
cbf_RBANS_data = merge( RBANS_data_BL, MD_data, by="CandID") 
plot(asl_avg_GM ~ total_scale_index_score_60_69, cbf_RBANS_data, xlab="RBANS total score", ylab="Grey matter CBF (ml/100g/min)")
abline(lm(asl_avg_GM ~ total_scale_index_score_60_69, cbf_RBANS_data))
summary(lm(asl_avg_GM ~ total_scale_index_score_60_69, cbf_RBANS_data))
plot(asl_PCN ~ total_scale_index_score_60_69, cbf_RBANS_data, xlab="RBANS total score", ylab="CBF in Precuneus (ml/100g/min)")
abline(lm(asl_PCN ~ total_scale_index_score_60_69, cbf_RBANS_data))
summary(lm(asl_PCN ~ total_scale_index_score_60_69, cbf_RBANS_data))


# Create super dataframe with variables of interest
CBF_MD_dataframe 	= CreateSuperDataframe()
CBF_MD_dataframe$Candidate_Age_Years = CBF_MD_dataframe$Candidate_Age / 12

# Age vs Gender graph
plot(asl_avg_GM ~ Candidate_Age_Years, females, col = "red", xlim=c(55,85), ylim=c(35,85))
points(asl_avg_GM ~ Candidate_Age_Years, males, col = "blue")
abline(rlm(asl_avg_GM ~ Candidate_Age_Years, females), col = "red")
abline(rlm(asl_avg_GM ~ Candidate_Age_Years, males),   col = "blue")
age_gender_model = CBF_MD_dataframe$asl_avg_GM ~ CBF_MD_dataframe$Candidate_Age_Years + CBF_MD_dataframe$Gender
summary(rlm(age_gender_model))
pvalue_age    = 2*min(1-pnorm(-3.1283), pnorm(-3.1283))  # 0.001758206 **
pvalue_gender = 2*min(1-pnorm(-2.2299), pnorm(-2.2299))	 # 0.02575408  *

# RBANS total score
plot(asl_avg_GM ~ total_scale_index_score_60_69, CBF_MD_dataframe)
abline(rlm( asl_avg_GM ~ total_scale_index_score_60_69, CBF_MD_dataframe))
summary(rlm( asl_avg_GM ~ total_scale_index_score_60_69, CBF_MD_dataframe))

###################
# USEFUL COMMANDS #
###################
# to get the list of available objects
ls( all.names = TRUE )
ls( pattern = "data" )

# to get the list of unique candID in prec_pass_data
unique(prec_pass_data$CandID)

# merge CBF and medication data into a new dataframe
med_use_cbf_pass = merge(curr_med_use_data,prec_pass_data,by=c("PSCID","CandID"))
candidate_info_cbf_pass = merge(candidate_info,prec_pass_data, by="CandID")
# select only subsample from a dataframe
T1_cbf_pass[T1_cbf_pass$bil_prec_WeightedAvgCBF<35,c(1,2,21,54)]

###### REMOVE FAILED VISITS FROM SPREADSHEETS
instrument_list = ls()
instrument_list = instrument_list[!instrument_list %in% c('LoadAdditionalSpreadsheet', 'LoadDataRelease')]
for (i in 1:length(instrument_list)) {
	
	# Initialize variables
	variable_name = instrument_list[i]
	
	assign( variable_name,
			RemoveFailedVisits(get(variable_name))
		  )

}




prec_pass_data = subset (cbf_data, Precuneus_QC=="pass" & Visit_label %in% c("NAPBL00","PREBL00"))

# load medication data (the spreadsheet below has been appended the changes noticed at baseline)
#curr_med_use_data = read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/for analyses/DATA_RELEASE_1.0_spreadsheets_BLonly_passPreCBF/current_med_use_with_changes_at_baseline.csv",h=T,sep=",")
# remove empty columns (find out a more efficient way to do this automatically)
#curr_med_use_data$Window_Difference <- NULL
#curr_med_use_data$X3_anticoagulant_su <- NULL
#curr_med_use_data$X3_anticoagulant_su_other <- NULL
#curr_med_use_data$X3_anticoagulant_su_bin_status <- NULL
#curr_med_use_data$X3_anticoagulant_prn <- NULL
#curr_med_use_data$X3_anticoagulant_prn_other <- NULL
#curr_med_use_data$X3_anticoagulant_prn_bin_status <- NULL
#curr_med_use_data$X4_corticosteroids_su <- NULL
#curr_med_use_data$X4_corticosteroids_su_other <- NULL
#curr_med_use_data$X4_corticosteroids_su_bin_status <- NULL
#curr_med_use_data$X5_reductase_inhibitors_prn <- NULL
#curr_med_use_data$X5_reductase_inhibitors_prn_other <- NULL
#curr_med_use_data$X5_reductase_inhibitors_prn_bin_status <- NULL
#curr_med_use_data$X10_osteoporosis_prn <- NULL
#curr_med_use_data$X10_osteoporosis_prn_other <- NULL
#curr_med_use_data$X10_osteoporosis_prn_bin_status <- NULL
#curr_med_use_data$X11_thyroid_prn <- NULL
#curr_med_use_data$X11_thyroid_prn_other <- NULL
#curr_med_use_data$X11_thyroid_prn_bin_status <- NULL
#curr_med_use_data$X13_diabetes_prn <- NULL
#curr_med_use_data$X13_diabetes_prn_other <- NULL
#curr_med_use_data$X13_diabetes_prn_bin_status <- NULL
#curr_med_use_data$X16_antidepressants_prn <- NULL
#curr_med_use_data$X16_antidepressants_prn_other <- NULL
#curr_med_use_data$X16_antidepressants_prn_bin_status <- NULL
#curr_med_use_data$X17_alpha_a1_prn <- NULL
#curr_med_use_data$X17_alpha_a1_prn_other <- NULL
#curr_med_use_data$X17_alpha_a1_prn_bin_status <- NULL





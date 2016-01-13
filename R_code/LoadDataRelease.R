# Load Data (data release, cbf and MD spreadsheets)
LoadSpreadsheets <- function (data_release_path, cbf_path, MD_path) {

	###### LOAD DATA RELEASE SPREADSHEETS
	LoadDataRelease( data_release_path )

	###### LOAD IMAGING DATA
	# load CBF data
	LoadAdditionalSpreadsheet( 'cbf', cbf_path )
	CleanUpCBF("cbf")

	# load MD data
	LoadAdditionalSpreadsheet( 'MD_data', MD_path )
	MD_data$apoe4 <- NULL

	# CLEAN UP DATA RELEASE SPREADSHEETS
	DataReleaseCleanUp( data_release_path )

}

CreateSuperDataframe      <- function () {

	# Create super dataframe with variables of interest
	CBF_MD_dataframe = merge ( 
		MD_data[,c("CandID", "asl_avg_GM", "asl_PCN", "MD_pcn")],
		cohort_candidate_info_data[,c("CandID", "Gender")],
		by="CandID",
		all.x = TRUE
	)

	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,
	   	general_information_data[,c("CandID","X4_education", "level_of_education")],
	   	by = "CandID",
		all.x = TRUE
	)

	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,
	   	subset( baseline_data, 
	   			!(PSCID == "MTL0033" & Visit_label == "PREBL00"),
	   			select = c( "CandID", 	 "Candidate_Age", 
	   					  "bp_systolic", "bp_diastolic",  "X8_bp_normal",  
	   					  "X9_pulse", 	 "X9_pulse_regularity" ) 
	   		  ),
	    by = "CandID",
		all.x = TRUE
	)

	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,
	    genetics_data[,c( "CandID", 			"E4_allele_Bin",  "K_variant_bin",    "BDNF_copie_bin", 
	    				  "Intron_M_allele_no", "TLR4_allele_no", "ppp2r1A_copie_no", "ApoE_HMGR_bin",  
	    				  "Reference_ApoE", 	"Reference_BDNF", "Reference_TLR4",   "Reference_ppp2r1a" )],
	    by = "CandID",
		all.x = TRUE
	)	
		
	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,																
		subset( lab_data, 
				Visit_label %in% c("NAPBL00", "PREBL00") & !(PSCID == "MTL0033" & Visit_label == "PREBL00"),
				select = c( "CandID", 
						    "X1_hb_result", 
						    "X2_ht_result", 
						    "X3_wbc_result",
						    "X4_rbc_result", 
						    "X5_platelet_result",
						    "X6_sodium_result",
						    "X7_potassium_result",
						    "X8_chloride_result",
						    "X9_urea_result",
						    "X10_creatinine_result",
						    "X11_protein_tot_result",
						    "X12_bilirubin_tot_result",
						    "X13_alk_phosph_result",
						    "X14_alt_result",
						    "X15_ast_result",
						    "X16_albumin_result",
						    "X17_ggt_result",
						    "X18_hba1c_result",
						    "X19_tsh_result", 
						    "X20_pt_result",
						    "X21_inr_result",
						    "X22_ptt_result" )
			  ),
	    by = "CandID",
		all.x = TRUE
	)	
	
	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,																
		subset( smell_identification_data, 
				Visit_label %in% c("NAPBL00", "PREBL00") & !(PSCID == "MTL0033" & Visit_label == "PREBL00"),
				select = c( "CandID",
					        "test_language",
					        "smoke_now",
					        "smoke_ever",
					        "ccsit_bsit",
					        "bsit_a",
					        "bsit_a_new",
					        "bsit_b",
					        "total_score",
					        "diagnosis")
			  ),
	    by = "CandID",
		all.x = TRUE
	)	
	
	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,																   
		subset(LP_data, Visit_label %in% c("NAPLP00", "PRELP00"))[,c( "CandID",
																	  "bp_sitting",
																	  "pulse_sitting",
																	  "ELISA_tau",
																	  "ELISA_b_amyloid",
																	  "ELISA_ptau",
																	  "tau_over_abeta_ratio",
																	  "ELISA_NFL",
																	  "ELISA_date",
																	  "X6plex_ApoAI",
																	  "X6plex_ApoAII",
																	  "X6plex_ApoB",
																	  "X6plex_ApoCII",
																	  "X6plex_ApoCIII",
																	  "X6plex_ApoE",
																	  "X6plex_date" )],
	    by = "CandID",
		all.x = TRUE
	)	
	
	CBF_MD_dataframe = merge ( 
		CBF_MD_dataframe,
		subset( RBANS_data, 
				Visit_label %in% c("NAPBL00", "PREBL00") & !(PSCID == "MTL0033" & Visit_label == "PREBL00"),
				select = c( "CandID",
						    "test_language",
						    "test_version",
						    "total_scale_index_score_60_69",
						    "immediate_memory_index_score_60_69",
						    "visuospatial_constructional_index_score_60_69",
						    "language_index_score_60_69",
						    "attention_index_score_60_69",
						    "delayed_memory_index_score_60_69" )
			  ),
		suffixes = c("_SNIF", "_RBANS"),
	    by = "CandID",
	    all.x = TRUE
	)
	
	return(CBF_MD_dataframe)

}

# Load Data Release Instruments's CSV spreadsheets
LoadDataRelease 		  <- function (path_to_data_release) 	   {

	# Create a list of spreadsheets found in directory "path_to_data_release"
	spreadsheet_list = list.files(path_to_data_release)
	
	# Loop through list of spreadsheets
	for (i in 1:length(spreadsheet_list)) {
		
		# Initialize variables
		variable_name = gsub(".csv", "", spreadsheet_list[i])
		spreadsheet   = paste(path_to_data_release, spreadsheet_list[i],sep="/")
		
		# Create instrument's dataframe from spreadsheet
		print(paste("Loading", spreadsheet_list[i]))
		assign( variable_name, 
			    read.table(spreadsheet, h=T, sep=",")	# read instrument spreadsheet
			  )
		
		# Convert CandID column of instrument's dataframe to a factor, unless instrument is DataDictionary
		if (variable_name != "DataDictionary") {
			assign( variable_name, 
				    transform(get(variable_name), CandID=as.factor(CandID)), # transform CandID to a factor
				    envir = .GlobalEnv # Option envir of assign allows global availability of the instrument's dataframe
				  )
		}
	}
}	

# Load additional spreadsheet like CBF and MD
LoadAdditionalSpreadsheet <- function (variable_name, spreadsheet) {
	
	# Create dataframe from spreadsheet
	print(paste("Loading", spreadsheet))
	assign( variable_name, read.table(spreadsheet, h=T, sep=",") )
	
	# Convert CandID column of the dataframe to a factor
	assign( variable_name, 
		    transform(get(variable_name), CandID=as.factor(CandID)), # transform CandID to a factor
		    envir = .GlobalEnv # Option envir of assign allows global availability of the instrument's dataframe
		  )
	
}

DataReleaseCleanUp 		  <- function (path_to_data_release) 	   {
	
	# Create a list of spreadsheets found in directory "path_to_data_release"
	spreadsheet_list = list.files(path_to_data_release)
	
	# Loop through list of spreadsheets
	for (i in 1:length(spreadsheet_list)) {
		# Initialize variables
		variable_name = gsub(".csv", "", spreadsheet_list[i]) 
		 
		if (variable_name == "AD8") 	 			 	 {	CleanUpAD8(variable_name) 			 } 
		if (variable_name == "baseline") 			 	 {	CleanUpBaseline(variable_name)		 }
		if (variable_name == "candidate_info") 		 	 {	CleanUpCandidateInfo(variable_name)  }
		if (variable_name == "cdr") 				 	 {	CleanUpCDR(variable_name)			 }
		if (variable_name == "enrollment") 			 	 {	CleanUpEnrollment(variable_name)	 }
		if (variable_name == "general_information")  	 {	CleanUpGeneralInfo(variable_name)	 }
		if (variable_name == "general_medical_history")  {	CleanUpGeneralMedHist(variable_name) }
		if (variable_name == "general_physical")		 { 	CleanUpGeneralPhys(variable_name)	 }
		if (variable_name == "genetics")				 {	CleanUpGenetics(variable_name)		 }
		if (variable_name == "handedness") 			 	 {	CleanUpHandedness(variable_name)	 }
	   	if (variable_name == "lab")						 {	CleanUpLab(variable_name)			 }
		if (variable_name == "LP")						 { 	CleanUpLP(variable_name)			 }
		if (variable_name == "MOCA")					 { 	CleanUpMOCA(variable_name)			 }
		if (variable_name == "RBANS")					 {	CleanUpRBANS(variable_name)			 }
		if (variable_name == "mri_feedbacks_adniT1") 	 {	CleanUpMRIfeedbacksT1(variable_name) }
		if (variable_name == "smell_identification")	 {	CleanUpOlfaction(variable_name)		 }
#		if (variable_name == "telephone_visit")			 {	CleanUpTelVisit(variable_name)		 }			
#		if (variable_name == "follow_up_visits") 	 	 {	CleanUpFollowUp(variable_name)		 }
# 		if (variable_name == "treatment_interruption") 	 {	CleanUpTreatmentInterr(variable_name)}
		if (variable_name == "mri_feedbacks_FinalnoRegQCedDTI") {	CleanUpMRIfeedbackDWI(variable_name) }
		
		
		# Remove failed visits
		exception_instruments = c( "DataDictionary",      "candidate_info", 
								   "inclusion_exclusion", "current_med_use", 
								   "adverse_events",      "telephone_visit", 
								   "follow_up_visits",    "treatment_interruption",
								   "mri_feedbacks_ASL",	  "mri_parameter_form",
								   "radiology_review",	  "ParticipantStatusHistory",
								   "tsi"
								 )
		if (variable_name %in% exception_instruments ) { 
			next 
		} else {
			new_variable_name = paste(variable_name,"data",sep="_")
			assign( variable_name,	RemoveFailedVisits(get(new_variable_name)) )
		}	
	}
	
}


CleanUpCBF 				  <- function (variable_name) {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "cbf_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID",                   "CandID",                 
 				 "SubprojectID",            "Visit_label",            
				 "Visit",					"MRI_Visit",                      
				 "CBF_QC",                  
				 "CBF_Comments",			"GM_AveragedCBF",			
				 "GM_QC",					"GM_Comments",            
				 "Precuneus_QC",            "Precuneus_Caveat",       
				 "Precuneus_Comments",      "bil_prec_WeightedAvgCBF", 
				 "bil_precl_StdDevCBF", 	"bil_prec_MinCBF",
				 "bil_prec_MaxCBF",			"bil_prec_VoxelCount",     
				 "bil_precl_VolSize",		"l_prec_WeightedAvgCBF",  
				 "l_prec_StdDevCBF",        "r_prec_WeightedAvgCBF",   
				 "r_prec_StdDevCBF"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}


CleanUpAD8 		  		  <- function (variable_name) 		       {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "AD8_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit_label", 
				 "Visit", 
				 "X1_judgment_problems", 
				 "X2_less_interest", 
				 "X3_repeat", 
				 "X4_trouble_learning_tool", 
				 "X5_forget_month_year", 
				 "X6_trouble_financial_affairs", 
				 "X7_trouble_appointments", 
				 "X8_daily_memory_trouble", 
				 "total_score"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpBaseline 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "baseline_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit_label",
				 "Visit",
				 "Date_taken",
				 "Candidate_Age",
				 "X5_glucose",
				 "X5_protein",
				 "bp_systolic",
				 "bp_diastolic",
				 "X8_bp_normal",
				 "X9_pulse",
				 "X9_pulse_regularity"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpCandidateInfo 	  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	cohort_name = "cohort_candidate_info_data"
	naproxen_name = "nap_candidate_info_data"
	
	# create two new dataframes (one for the cohort, one for the naproxen group)
	assign( cohort_name,
			subset( get(variable_name), SubprojectID == "PreventAD" ),
			envir = .GlobalEnv
		  )
	assign( naproxen_name,
		  	subset( get(variable_name), SubprojectID == "NaproxenTrial"),
		  	envir = .GlobalEnv
		  )
	
}

CleanUpCDR 				  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "cdr_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID",
				 "Visit_label", 
				 "Visit",
				 "cdr"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpEnrollment 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "enrollment_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID",
				 "Visit_label", 
				 "Visit",
				 "Date_taken",
				 "Candidate_Age",
				 "X1_hospitalization",
				 "X2_new_symptoms",
				 "X2_new_conditions",
				 "X2_worsened_condition",
				 "X3_medication_changes",
				 "X6_EKG_result",
				 "X6_EKG_abnormal",
				 "bp_systolic",
				 "bp_diastolic",
				 "X8_bp_normal",
				 "X9_pulse",
				 "X9_pulse_regularity"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}
 
CleanUpGeneralInfo 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name   = "general_information_data"
	family_history_name = "family_history_data" 
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID",
				 "Visit_label", 
				 "Visit",
				 "Date_taken",
				 "Candidate_Age",
				 "X2_ethnicity",
				 "ethnicity_other",
				 "X3_marital_stat",
				 "X5_work",
				 "X6_retired",
				 "X7_language",
				 "X8_residence",
				 "X9_living_situation",
				 "X11_appetite_habits",
				 "abnormal_appetite",
				 "X12_sleeping_habits",
				 "abnormal_sleeping",
				 "X13_cigarette_smoking",
				 "cigarette_stopped",
				 "X14_cigar_smoking",
				 "cigar_stopped",
				 "X15_pipe_smoking",
				 "pipe_stopped",
				 "X16_drug_consumption",
				 "drug_stopped",
				 "X17_alcohol_consumption",
				 "alcohol_stopped",
				 "X18_dependency",
				 "X21_hobbies",
				 "X4_education",
				 "level_of_education",
				 "X22_activity_level" 
			   )
	
	family_hist_fields = c( "PSCID", 
				 			"CandID",
				 			"Visit_label", 
				 			"Visit",
				 			"X10_family_history_alzheimer",
				 			"family_member_1",
				 			"family_member_death_age_1",
				 			"family_member_age_onset_1",
							"family_member_2",
				 			"family_member_death_age_2",
				 			"family_member_age_onset_2",
				 			"family_member_3",
				 			"family_member_death_age_3",
				 			"family_member_age_onset_3",
				 			"family_member_4",
				 			"family_member_death_age_4",
				 			"family_member_age_onset_4",
				 			"family_member_5",
				 			"family_member_death_age_5",
				 			"family_member_age_onset_5",
				 			"family_member_6",
				 			"family_member_death_age_6",
				 			"family_member_age_onset_6",
				 			"family_member_7",
				 			"family_member_death_age_7",
				 			"family_member_age_onset_7",
				 			"X10a_number_siblings",
				 			"X10b_other_family_dementia_history"
						  )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	assign( family_history_name,
			get(variable_name)[,family_hist_fields],
			envir = .GlobalEnv
		  )
	
	# create binary for family types
	family_type = c("father","mother", "sibling")
	for (i in 1:length(family_type)) {
		family_member = family_type[i]
		new_field_name = paste(family_member,"bin",sep="_")	
		assign( new_field_name,
				ifelse( (!is.na(family_history_data$family_member_1) & family_history_data$family_member_1 == family_member) 
					  | (!is.na(family_history_data$family_member_2) & family_history_data$family_member_2 == family_member) 
				  	  | (!is.na(family_history_data$family_member_3) & family_history_data$family_member_3 == family_member)
					  | (!is.na(family_history_data$family_member_4) & family_history_data$family_member_4 == family_member)
					  | (!is.na(family_history_data$family_member_5) & family_history_data$family_member_5 == family_member)
					  | (!is.na(family_history_data$family_member_6) & family_history_data$family_member_6 == family_member)
					  | (!is.na(family_history_data$family_member_7) & family_history_data$family_member_7 == family_member),
					  	"yes",
					  	"no"
					  )
			  )
		assign( new_field_name, 
				as.data.frame(get(new_field_name)), 
				envir = .GlobalEnv
			  )
	}
	names(father_bin)  = "father_bin"
	names(mother_bin)  = "mother_bin"
	names(sibling_bin) = "sibling_bin"
	assign( family_history_name, cbind(get(family_history_name), father_bin ) )
	assign( family_history_name, cbind(get(family_history_name), mother_bin ) )
	assign( family_history_name, cbind(get(family_history_name), sibling_bin) )
	# create a parent type
	parent_bin = ifelse( family_history_data$father_bin == "yes" | family_history_data$mother_bin == "yes", "yes", "no")
	assign( family_history_name, cbind(get(family_history_name), parent_bin ) )
	
	# compute average age of family onset
	family_history_data$family_member_age_onset_1 = as.character(family_history_data$family_member_age_onset_1)
	family_history_data$family_member_age_onset_2 = as.character(family_history_data$family_member_age_onset_2)
	family_history_data$family_member_age_onset_3 = as.character(family_history_data$family_member_age_onset_3)
	family_history_data$family_member_age_onset_4 = as.character(family_history_data$family_member_age_onset_4)
	family_history_data$family_member_age_onset_5 = as.character(family_history_data$family_member_age_onset_5)
	family_history_data$family_member_age_onset_6 = as.character(family_history_data$family_member_age_onset_6)
	family_history_data$family_member_age_onset_7 = as.character(family_history_data$family_member_age_onset_7)
	family_history_data$family_member_age_onset_1 = as.numeric(family_history_data$family_member_age_onset_1)
	family_history_data$family_member_age_onset_2 = as.numeric(family_history_data$family_member_age_onset_2)
	family_history_data$family_member_age_onset_3 = as.numeric(family_history_data$family_member_age_onset_3)
	family_history_data$family_member_age_onset_4 = as.numeric(family_history_data$family_member_age_onset_4)
	family_history_data$family_member_age_onset_5 = as.numeric(family_history_data$family_member_age_onset_5)
	family_history_data$family_member_age_onset_6 = as.numeric(family_history_data$family_member_age_onset_6)
	family_history_data$family_member_age_onset_7 = as.numeric(family_history_data$family_member_age_onset_7)
	family_onset_matrix = family_history_data[,c("family_member_age_onset_1", "family_member_age_onset_2", "family_member_age_onset_3", "family_member_age_onset_4", "family_member_age_onset_5", "family_member_age_onset_6", "family_member_age_onset_7")]
	family_mean_onset = as.data.frame( rowMeans(family_onset_matrix, na.rm = TRUE) )
	names( family_mean_onset ) = "family_mean_age_onset"

	# compute number of family member with AD
	family_number_matrix = family_history_data[,c("family_member_1", "family_member_2", "family_member_3", "family_member_4", "family_member_5", "family_member_6", "family_member_7")]
	number_of_family_with_AD = as.data.frame( apply(family_number_matrix, 1, function(x) length(which(!is.na(x)))) )
	names( number_of_family_with_AD ) = "number_of_family_member_with_AD"
	assign( family_history_name, cbind(get(family_history_name), family_mean_onset, number_of_family_with_AD), envir = .GlobalEnv )

}

CleanUpGeneralMedHist 	  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "general_medical_history_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID",
				 "Visit_label", 
				 "Visit",
				 "X1_peptic_ulcer_complications",               
				 "X2_peptic_ulcer_uncomplicated",
				 "X3_treatment_hypertension",                   
				 "X4_treatment_hyperlipidemia",
				 "X5_treatment_diabetes",                       
				 "X6_strokes",
				 "X6_strokes_number",                           
				 "X6_strokes_year",                             
				 "X7_tia",                                      
				 "X8_head_injury",
				 "X9_past_conditions",                          
				 "X10_past_diagnoses",
				 "X10_a_cancer_specify",                        
				 "X10_b_depression_specify",                    
				 "X10_e_myocardial_infraction_specify",         
				 "X10_j_other_neurodegenerative_specify",       
				 "X10_k_epilepsy_specify",                      
				 "X10_q_autoimmune_disorder_specify",           
				 "X10_v_other_specify",                         
				 "X11_prostate_hypertrophy",                    
				 "X12_other_diagnoses",
				 "X12_d_other_liver_specify",                   
				 "X13_other_diagnoses",                  
				 "X13_c_other_kidney_specify",
				 "X14_other_diagnoses",
				 "X15_surgery",                           
				 "X15_surgery_1",                     
				 "X15_surgery_1_year",              
				 "X16_hospitalizations",
				 "X16_hospitalizations_number"                
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
		
}

CleanUpGeneralPhys 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "general_physical_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID",
				 "Visit_label", 
				 "Visit",
				 "X1_blood_pressure",                             
				 "X2_blood_pressure_range",
 			     "X3_pulse",   
 			     "X3_pulse_regularity",
 			     "X4_weight",
 			     "X5_height",
 			     "bmi",
 			     "X7_glucose",
 			     "X7_protein",
 			     "X9_bruising",
 			     "X9_peripheral_oedema",
 			     "X9_lower_back_pain",
 			     "X9_tinnitus",
 			     "X9_dizziness_vertigo",
 			     "X9_insomnia",
 			     "X9_headache",
 			     "X9_numbness",
 			     "X9_tremor",
 			     "X10_general_appearance",
 			     "X10_skin",
 			     "X10_eyes",
 			     "X10_ears",
 			     "X10_mouth",
 			     "X10_heart",
 			     "X10_lung",
 			     "X10_abdomen",
 			     "X10_musculoskeletal",
 			     "X13_exclusionary_conditions",
 			     "cranial_II",
 			     "cranial_perral",
 			     "cranial_visual_fields",
 			     "cranial_III_IV_VI",
 			     "cranial_EOM",
 			     "cranial_nystagmus",
 			     "cranial_V",
 			     "cranial_VII",
 			     "cranial_VIII",
 			     "cranial_IX_X",
 			     "cranial_XI",
 			     "cranial_XII",
 			     "motor_arms_proximal",
 			     "motor_arms_distal",
 			     "motor_legs_proximal",
 			     "motor_legs_distal",
 			     "motor_tone",
 			     "motor_bulk",
 			     "sensory_pinprick",
 			     "sensory_fine_touch",
 			     "sensory_position",
 			     "sensory_vibration",
 			     "cerebellar_finger_nose",
 			     "cerebellar_heel_shin",
 			     "frontal_lobe_suck",
 			     "frontal_lobe_snout",
 			     "frontal_lobe_grasp",
 			     "frontal_lobe_palmo_mental",
 			     "frontal_lobe_glabellar",
 			     "gait_standard",
 			     "gait_tandem",
 			     "gait_rhomberg",
 			     "biceps_right",
 			     "biceps_left",
 			     "triiceps_left",
 			     "triceps_right",
 			     "brachioradialis_right",
 			     "brachioradialis_left",
 			     "knee_right",
 			     "knee_left",
 			     "ankle_right",
 			     "ankle_left",
 			     "plantar_response_right",
 			     "plantar_response_left",
 			     "X14_cranial_nerve",
 			     "X14_tendon_reflexes",
 			     "X14_motor_exam",
 			     "X14_gait_station",
 			     "X14_primitive_reflexes",
 			     "X14_sensorium",
 			     "X14_range_affect",
 			     "X14_expressive_speech"        
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpGenetics 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "genetics_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "ApoE",
				 "apoE_allele_no",
				 "E4_allele_Bin",
				 "Reference_ApoE",
				 "BchE_K_variant",
				 "K_variant_copie_no",
				 "K_variant_bin",
				 "BDNF",
				 "BDNF_allele_no",
				 "BDNF_copie_bin",
				 "Reference_BDNF",
				 "HMGR_Intron_M",
				 "Intron_M_allele_no",
				 "TLR4_rs_4986790",
				 "TLR4_allele_no",
				 "Reference_TLR4",
				 "PPP2r1A_rs_10406151",
				 "ppp2r1A_allele_no",
				 "ppp2r1A_copie_no",
				 "Reference_ppp2r1a"
			   )
	to_factorize = c()
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
	# factorize genetics data
	genetics_data$apoE_allele_no 	 = as.factor(genetics_data$apoE_allele_no	 )
	genetics_data$E4_allele_Bin 	 = as.factor(genetics_data$E4_allele_Bin	 )
	genetics_data$K_variant_copie_no = as.factor(genetics_data$K_variant_copie_no)
	genetics_data$K_variant_bin 	 = as.factor(genetics_data$K_variant_bin	 )
	genetics_data$BDNF_allele_no 	 = as.factor(genetics_data$BDNF_allele_no 	 )
	genetics_data$BDNF_copie_bin 	 = as.factor(genetics_data$BDNF_copie_bin	 )
	genetics_data$Intron_M_allele_no = as.factor(genetics_data$Intron_M_allele_no)
	genetics_data$TLR4_allele_no 	 = as.factor(genetics_data$TLR4_allele_no	 )
	genetics_data$ppp2r1A_allele_no  = as.factor(genetics_data$ppp2r1A_allele_no )
	genetics_data$ppp2r1A_copie_no   = as.factor(genetics_data$ppp2r1A_copie_no  )
	# HMGR is protective and cancel ApoE4 effect.
	genetics_data$ApoE_HMGR_bin = as.factor( ifelse( genetics_data$E4_allele_Bin == 1 & genetics_data$Intron_M_allele_no != 2, 
											 		 1, 
											 		 0 
									 	   		   ) 
							 			   )
	assign( new_variable_name,
			get(new_variable_name),
			envir = .GlobalEnv
		  )

}

CleanUpHandedness 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "handedness_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit",
				 "result",
				 "interpretation"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpLab 				  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "lab_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID",				 		"CandID", 
				 "Visit_label",						"Visit",					
				 "Date_taken",					"Candidate_Age",
				 "X1_hb_result",   			    "X2_ht_result",
				 "X3_wbc_result",			    "X4_rbc_result",
				 "X5_platelet_result",			"X6_sodium_result",
				 "X7_potassium_result",		    "X8_chloride_result",
				 "X9_urea_result",			    "X10_creatinine_result",
				 "X11_protein_tot_result",	 	"X12_bilirubin_tot_result",
				 "X13_alk_phosph_result",	    "X14_alt_result",
				 "X15_ast_result",			 	"X16_albumin_result",
				 "X17_ggt_result",			 	"X18_hba1c_result",
				 "X19_tsh_result",				"X20_pt_result",
				 "X21_inr_result",				"X22_ptt_result"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	lab_data_names = names(lab_data)	  
	
	# calculate whether in normal range
	list_test = c( "X1_hb",				   "X2_ht",
				   "X3_wbc",			   "X4_rbc",
				   "X5_platelet",		   "X6_sodium",
				   "X7_potassium",		   "X8_chloride",
				   "X9_urea",			   "X10_creatinine",
				   "X11_protein_tot",	   "X12_bilirubin_tot",
				   "X13_alk_phosph",	   "X14_alt",
				   "X15_ast",			   "X16_albumin",
				   "X17_ggt",			   "X18_hba1c",
				   "X19_tsh",			   "X20_pt",
				   "X21_inr",			   "X22_ptt"
				 )
				 
	for (i in 1:length(list_test)) {
		# 
		tmp = ComputeLabInterpretation(list_test[i], "lab")
		
		assign( new_variable_name, cbind(lab_data, tmp[,5]) )
		interpretation_name	= paste(list_test[i], "interpretation", sep="_")
		lab_data_names    = c(lab_data_names, interpretation_name)
		names( lab_data ) = lab_data_names
		assign( new_variable_name, get(new_variable_name), envir = .GlobalEnv )
	}		  

}

CleanUpLP 				  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "LP_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID",				 			"CandID", 
				 "Visit_label",						"Visit",					
				 "Date_taken",						"Candidate_Age",
				 "nicotine_drugs",					"new_medications_dosages",         
				 "illnesses_surgeries_accidents",	"proceed", 
 				 "bp_sitting",                  		"pulse_sitting",
 				 "X1_fluid_result",					"X2_rbc_result",
 				 "X3_neutrophils_result",           "X4_lymphocytes_result",           
 				 "X5_mononuclears_result", 			"X6_epithelial_cells_result",			
 				 "X8_wbc_result",                                 
 				 "ELISA_tau",                       "ELISA_b_amyloid",                 
 				 "ELISA_ptau",                      "ELISA_date",                       
				 "ELISA_NFL",                       "X6plex_ApoAI",                     
				 "X6plex_ApoAII",                   "X6plex_ApoB",                      
				 "X6plex_ApoCII",                   "X6plex_ApoCIII",                   
				 "X6plex_ApoE",                     "X6plex_date",
				 "X9_glucose_result",				"X10_micro_protein_result"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	LP_data_names = names(LP_data)	  

	# calculate whether in normal range
	list_test = c( "X9_glucose",	"X10_micro_protein" )
	
	for (i in 1:length(list_test)) {
		# 
		tmp = ComputeLabInterpretation(list_test[i], "LP")
		
		assign( new_variable_name, cbind(LP_data, tmp[,5]) )
		interpretation_name	= paste(list_test[i], "interpretation", sep="_")
		LP_data_names = c(LP_data_names, interpretation_name)
		names( LP_data ) = LP_data_names
		assign( new_variable_name, get(new_variable_name), envir = .GlobalEnv )
	}
	
	LP_data$tau_over_abeta_ratio = LP_data$ELISA_tau/LP_data$ELISA_b_amyloid
	assign( new_variable_name, get(new_variable_name), envir = .GlobalEnv )
	
}

CleanUpMRIfeedbacksT1 	  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "mri_feedbacks_adniT1_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit",
				 "MRI_Visit",
				 "Susceptibility_Artefact_Dental_Work",
				 "Susceptibility_Artefact_Anatomy",
				 "Susceptibility_Artefact_Ear_Canals"			   
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpMRIfeedbackDWI	  <- function (variable_name)			   {
	
		# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "mri_feedbacks_FinalnoRegQCedDTI_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit",
				 "MRI_Visit",
				 "File_QC_status",
				 "File_Caveat",
				 "Color_Artifact",
				 "Red_Artifact",
				 "Intensity_Comment",
				 "Movement",
				 "Gradient_Wise_Motion",
				 "Motion_Comment",
				 "Medium_AP_Wrap",
				 "Top_Brain_Cut_Off",
				 "Too_Few_Remaining_Gradients"			   
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
	
}

CleanUpMOCA 			  <- function (variable_name) 			   { 
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "MOCA_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit",
				 "visuospatial_score",
				 "naming_score",
				 "attention_numbers_score",
				 "attention_letters_score",
				 "attention_subtractions_score",
				 "attention_total_score",
				 "language_repeat_score",
				 "language_fluency_score",
				 "language_total_score",
				 "abstraction_score",
				 "delayed_recall_score",
				 "orientation_score",
				 "total_score"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )

}

CleanUpRBANS 			  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "RBANS_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 									  			  "CandID", 
				 "Visit_label",							    			  "Visit",
				 "test_language",										  "test_version",
				 "immediate_memory_index_score_60_69",				 	  "visuospatial_constructional_index_score_60_69",
	 			 "language_index_score_60_69",                      	  "attention_index_score_60_69",                          
				 "delayed_memory_index_score_60_69",				 	  "total_scale_index_score_60_69",                        
				 "immediate_memory_index_score",                          "visuospatial_constructional_index_score",              
				 "language_index_score",                                  "attention_index_score",                                
				 "delayed_memory_index_score",                            "total_scale_index_score",
				 "list_learning_primacy_trial1",                          "list_learning_mid_trial1",                             
				 "list_learning_recency_trial1",                          "list_learning_primacy",                                
				 "list_learning_mid",                                     "list_learning_recency",                                
				 "list_learning_repetitions",                             "list_learning_intrusions", 
				 "story_memory_primacy_trial1",							  "story_memory_mid_trial1",
				 "story_memory_recency_trial1",							  "story_memory_primacy",
				 "story_memory_mid",									  "story_memory_recency",
				 "story_memory_intrusions",								  "story_memory_repetitions",
				 "semantic_fluency_intrusions",                           "semantic_fluency_repetitions",
				 "list_recall_intrusions",								  "list_recall_repetitions",                               
				 "list_recall_primacy",									  "list_recall_mid",                                       
				 "list_recall_recency",									  "list_recall_failure_2",                                 
				 "list_recall_failure_3",								  "list_recall_failure_4",                                
			     "list_recognition20",									  "list_recognition_score",                               
			     "list_recognition_targets_score",						  "list_recognition_distractors_score",                    		
			     "list_recognition_failure_2",							  "list_recognition_failure_3",                            
			     "list_recognition_failure_4",
			     "story_recall_score",                                    "story_recall_primacy",                                 
				 "story_recall_mid",                                      "story_recall_recency",                                
				 "story_recall_intrusions",                               "story_recall_repetitions" ,                            
				 "story_recall_failure_2",                                "list_recognition19",                                   
				 "verbal_memory_total_intrusions"  ,                      "verbal_memory_total_repetitions" 
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

CleanUpOlfaction 		  <- function (variable_name) 			   {
	
	# create variable with name of new dataframe that will be created (AD8_scores)
	new_variable_name = "smell_identification_data"
	
	# create variable with fields to keep
	to_keep = c( "PSCID", 
				 "CandID", 
				 "Visit_label",
				 "Visit",
				 "test_done",
				 "test_language",
				 "perceive_weak_normal",
				 "smoke_now",
				 "smoke_now_cigarettes_day",
				 "smoke_now_cigars_day",
				 "smoke_now_other",
				 "smoke_inhale",
				 "smoke_now_start",
				 "smoke_ever",
				 "smoke_ever_start",
				 "smoke_ever_stop",
				 "smoke_ever_cigarettes_day",
				 "smoke_ever_cigars_day",
				 "smoke_ever_other",
				 "smoke_smell_ability",
				 "ccsit_bsit",
				 "bsit_a",
				 "bsit_a_new",
				 "bsit_b",
				 "total_score",
				 "diagnosis"
			   )
	
	# create new dataframe with relevant fields from to_keep
	assign( new_variable_name,
			get(variable_name)[,to_keep],
			envir = .GlobalEnv
		  )
	
}

ComputeLabInterpretation  <- function (test_name, instrument_name) {
	
		# initialize variables
		range_variable_name = paste(test_name, "normal_range",   sep="_")
		value_variable_name = paste(test_name, "result", 	     sep="_")
		min_variable_name   = paste(test_name, "min", 	   	     sep="_")
		max_variable_name   = paste(test_name, "max", 	   	     sep="_")
		interpretation_name	= paste(test_name, "interpretation", sep="_")
		
		print(paste("Creating",interpretation_name))
		
		# Create tmp dataframe with ranges and values
		tmp 		  = get(instrument_name)[,c(value_variable_name, range_variable_name)]
		tmp[,2] 	  = as.character(tmp[,2])
		tmp[,1]		  = as.character(tmp[,1])
		tmp$min 	  = NA
		tmp$max 	  = NA
		tmp$interpret = NA
		names(tmp) = c(value_variable_name, range_variable_name, min_variable_name, max_variable_name, interpretation_name) 
		
		# Split ranges column into min and max values and compute the interpretation for each line
		for (line in 1:length(tmp[,2])) { 
			tmp[line,c(3,4)] = unlist(strsplit(tmp[line,2],"-"))
			tmp[line,5] = ifelse( as.numeric(tmp[line,1]) > as.numeric(tmp[line,3]) & as.numeric(tmp[line,1]) < as.numeric(tmp[line,4]) ,"normal","abnormal")
		}
		
		# return the tmp dataframe
		return(tmp)
	
}


# Remove Failed visits from dataframes
RemoveFailedVisits 		  <- function (dataframe) 				   {
	
	# removing Visit and MRI_Visit that are set to Fail from dataframe
	if ("Visit" %in% names(dataframe))     { dataframe = subset(dataframe, !Visit %in% c("Fail","Failure"))	}
	if ("MRI_Visit" %in% names(dataframe)) { dataframe = subset(dataframe, MRI_Visit != "Fail") }
	
	return(dataframe)
	
}



#CleanUpFollowUp <- function(variable_name) {
#	
#	# create variable with name of new dataframe that will be created (AD8_scores)
#	new_variable_name = "follow_up_visits_data"
#	
#	# create variable with fields to keep
#	to_keep = c( "PSCID", 
#				 "CandID",
#				 "Visit_label", 
#				 "Visit",
#				 "Date_taken",
#				 "Candidate_Age",
#				 "X1_hospitalization",
#				 "X2_medication_changes",
#				 "X8_EKG_result",
#				 "X8_EKG_abnormal",
#				 "X9_systolic_bp",
#				 "X9_diastolic_bp",
#				 "X10_bp_range",
#				 "X11_pulse",
#				 "X11_pulse_regularity",
#				 "X12_weight",
#				 "X12_weight_units",			
#			   )
#	
#	# create new dataframe with relevant fields from to_keep
#	assign( new_variable_name,
#			get(variable_name)[,to_keep],
#			envir = .GlobalEnv
#		  )
#	
#}

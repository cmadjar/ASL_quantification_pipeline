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

### look at data that are passed for both MD and CBF data
cbf_general_physical = merge( general_physical_data, MD_data, by="CandID")
 
hist









###################
# USEFUL COMMANDS #
###################
# to get the list of available objects
ls(all.names = TRUE)

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





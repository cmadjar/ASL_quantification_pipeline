
# clear workspace
#rm(list = ls())

# Load visual QC spreadsheet
#CBF_visualQC = read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/QC spreadsheets/DKT_template/visual_QC_spreadsheet_20150928.csv", h=T, sep=",")

# Load CBF averages spreadsheet
CBF_avg = read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/QC spreadsheets/DKT_template/GM_avg_and_DKT_CBF_2015-11-28_without_failed_visual_QC.csv", h=T, sep=",")

# Load mincstats spreadsheet
mncstat = read.table("~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/QC spreadsheets/DKT_template/GM_avg_and_DKT_CBF_2015-11-30_minstats.csv", h=T, sep=",")

# count number of NAs per row (aka subject)
apply(CBF_avg, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

# Run automated QC outlier technique
#automated_qc_outliers_technique(CBF_avg_QC)

### This is not working. Need to find out why when have time to do so
remove_failed_visual_QC_visits <- function(CBF_visualQC, CBF_avg) {
	
	# get the failed visual QC data only (when either CBF map failed or registration failures)
	failed_visual_QC = subset(CBF_visualQC[,c(2,4,11,15)], GM_QC == "fail" | CBF_QC == "fail")
#	write.csv(failed_visual_QC, file="~/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/QC/QC spreadsheets/DKT_template/failed_visual_QC_CBF.csv")

	# replace all values in row of subjects that failed visual QC
	for (i in 1:dim(failed_visual_QC)[1]) {
		# grep the candidate and visit
		candID = failed_visual_QC[i,1]
		visit  = failed_visual_QC[i,2]
		
		# replace all extracted CBF by NA for that candidate/visit
		CBF_avg[ CBF_avg$CandID == candID & CBF_avg$Visit_label == visit, 3:dim(CBF_avg)[2] ] = NA
	}
	
	return(CBF_avg)
	
}


# Get list of ROIs
automated_qc_outliers_technique <- function(CBF_avg) {
	
	# Loop through all columns
	for (col_nb in 1:dim(CBF_avg)[2]) {
		# Go to the next column if column name is not CBF average
		if (names(CBF_avg[col_nb]) == "CandID" 
			|| names(CBF_avg[col_nb]) == "Visit" 
			|| names(CBF_avg[col_nb]) == "CBF_map"
		   ) {
			next
		} else {
			# compute lower and upper limits to use to find outliers for that columns
			lower_CBF_limit = mean(CBF_avg[, col_nb], na.rm = T) - sd(CBF_avg[, col_nb], na.rm = T)
			upper_CBF_limit = mean(CBF_avg[, col_nb], na.rm = T) + sd(CBF_avg[, col_nb], na.rm = T)
			# find the outlier values and replace them by NA
			CBF_avg = find_outliers(CBF_avg, lower_CBF_limit, upper_CBF_limit, col_nb)
		}
		
	}
	
	return(CBF_avg)
	
}



# Replace outliers by NA
find_outliers <- function(CBF_avg, lower_limit, upper_limit, col_nb) {
	
	# loop through rows
	for (row_nb in 1:dim(CBF_avg)[1]) {
		# replace outliers by NA
		if (is.na(CBF_avg[row_nb,col_nb])) {
			next
		} else if (CBF_avg[row_nb, col_nb] >= upper_limit) {
			as.character(CBF_avg[row_nb, col_nb])
			CBF_avg[row_nb, col_nb] = NA
		} else if (CBF_avg[row_nb, col_nb] <= lower_limit) {
			as.character(CBF_avg[row_nb, col_nb])
			CBF_avg[row_nb, col_nb] = NA
		} 
	}

	return(CBF_avg)	
		
}
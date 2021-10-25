## ERCC normalization of ERCC Spike-In RNAs

## Import library (tidyverse)
suppressMessages(library(tidyverse))

## get inputs <UTS output> <uniquely mapped reads> <dilution factor>  <ERCC Mix>
args = commandArgs(trailingOnly = TRUE)

## Set variables
path_abundance = args[1]
mapped_reads = as.numeric(args[2])
dilution = as.numeric(args[3])
ERCC = as.numeric(args[4])

## Get sample ID
sample_id = gsub(".genes.results", "", basename(path_abundance))

## Output file path and name
output_path = file.path(dirname(path_abundance), paste0("ERCCnorm_", sample_id, ".genes.results"))

## Import counts file
abundance = as.data.frame(suppressMessages(read_table(path_abundance)))

## Make a matrix using only unfiltered ERCC genes
ERCC_abundance = abundance[grepl("ERCC-", abundance[[1]]),]

## Make a matrix using only unfiltered non-ERCC genes
abundance = abundance[!grepl("ERCC-", abundance[[1]]),]

## Note: all samples should contain ERCC spike-in RNAs and thus ERCC counts.
exit <- function() { invokeRestart("abort") } # set exit function
if (colSums(ERCC_abundance[4]) == 0) {
	cat("Sample", sample_id, "had no detectable ERCC spike-ins\n")
	exit()
}

## Calculate log-base-2 of the ERCC RPKM and non-ERCC RPKM
ERCC_abundance$log2RPKM <- log2(ERCC_abundance[[5]] + 1)
abundance$log2RPKM <- log2(abundance[[5]] + 1)

##### Create a standard curve that relates ERCC Spike-In RNA-seq RPKM to
##### absolute amount of each ERCC transcript added

## Prepare ERCC concentrations
ERCC_Mix = as.data.frame(1:92)

## Get ERCC gene IDs
ERCC_Mix$gene_id <- c("ERCC-00130", "ERCC-00004", "ERCC-00136", "ERCC-00108", "ERCC-00116", "ERCC-00092", "ERCC-00095", "ERCC-00131", 
		         "ERCC-00062", "ERCC-00019", "ERCC-00144", "ERCC-00170", "ERCC-00154", "ERCC-00085", "ERCC-00028", "ERCC-00033",
		         "ERCC-00134", "ERCC-00147", "ERCC-00097", "ERCC-00156", "ERCC-00123", "ERCC-00017", "ERCC-00083", "ERCC-00096", 
		         "ERCC-00171", "ERCC-00009", "ERCC-00042", "ERCC-00060", "ERCC-00035", "ERCC-00025", "ERCC-00051", "ERCC-00053",
		         "ERCC-00148", "ERCC-00126", "ERCC-00034", "ERCC-00150", "ERCC-00067", "ERCC-00031", "ERCC-00109", "ERCC-00073",
	                 "ERCC-00158", "ERCC-00104", "ERCC-00142", "ERCC-00138", "ERCC-00117", "ERCC-00075", "ERCC-00074", "ERCC-00113", 
	                 "ERCC-00145", "ERCC-00111", "ERCC-00076", "ERCC-00044", "ERCC-00162", "ERCC-00071", "ERCC-00084", "ERCC-00099",
	                 "ERCC-00054", "ERCC-00157", "ERCC-00143", "ERCC-00039", "ERCC-00058", "ERCC-00120", "ERCC-00040", "ERCC-00164",
	                 "ERCC-00024", "ERCC-00016", "ERCC-00012", "ERCC-00098", "ERCC-00057", "ERCC-00002", "ERCC-00046", "ERCC-00003",
	                 "ERCC-00043", "ERCC-00022", "ERCC-00112", "ERCC-00165", "ERCC-00079", "ERCC-00078", "ERCC-00163", "ERCC-00059",
	                 "ERCC-00160", "ERCC-00014", "ERCC-00077", "ERCC-00069", "ERCC-00137", "ERCC-00013", "ERCC-00168", "ERCC-00041",
	                 "ERCC-00081", "ERCC-00086", "ERCC-00061", "ERCC-00048")

## Get ERCC spike-in RNA concentrations
if (ERCC == 1) {
	  ## Get concentrations of ERCC Spike-In RNAs Mix 1
	  ERCC_Mix$attomoles_per_ul <- c(30000, 7500, 1875, 937.5, 468.75, 234.375, 117.1875, 117.1875, 58.59375, 29.296875, 29.296875, 
					 14.6484375, 7.32421875, 7.32421875, 3.66210938, 1.83105469, 1.83105469, 0.91552734, 0.45776367,
					 0.45776367, 0.22888184, 0.11444092, 0.02861023, 15000, 3750, 937.5, 468.75, 234.375, 117.1875,
					 58.59375, 58.59375, 29.296875, 14.6484375, 14.6484375, 7.32421875, 3.66210938, 3.66210938,
					 1.83105469, 0.91552734, 0.91552734, 0.45776367, 0.22888184, 0.22888184, 0.11444092, 0.05722046, 
					 0.01430512, 15000, 3750, 937.5, 468.75, 234.375, 117.1875, 58.59375, 58.59375, 29.296875,
					 14.6484375, 14.6484375, 7.32421875, 3.66210938, 3.66210938, 1.83105469, 0.91552734, 0.91552734,
					 0.45776367, 0.22888184, 0.22888184, 0.11444092, 0.05722046, 0.01430512, 15000, 3750, 937.5, 468.75,
					 234.375, 117.1875, 58.59375, 58.59375, 29.296875, 14.6484375, 14.6484375, 7.32421875, 3.66210938,
					 3.66210938, 1.83105469, 0.91552734, 0.91552734, 0.45776367, 0.22888184, 0.22888184, 0.11444092,
					 0.05722046, 0.01430512)
	} else if (ERCC == 2) {
	  # Get concentrations of ERCC Spike-In RNAs Mix2
	  ERCC_Mix$attomoles_per_ul <- c(7500, 1875, 468.75, 234.375, 117.1875, 58.59375, 29.296875, 29.296875, 14.6484375, 7.3242187,
					 7.32421875, 3.66210938, 1.83105469, 1.83105469, 0.91552734, 0.45776367, 0.45776367, 0.22888184,
					 0.11444092, 0.11444092, 0.05722046, 0.02861023, 0.00715256, 15000, 3750, 937.5, 468.75, 234.375,
					 117.1875, 58.59375, 58.59375, 29.296875, 14.6484375, 14.6484375, 7.32421875, 3.66210938, 3.66210938,
					 1.83105469, 0.91552734, 0.91552734, 0.45776367, 0.22888184, 0.22888184, 0.11444092, 0.05722046, 
					 0.01430512, 22500, 5625, 1406.25, 703.125, 351.5625, 175.78125, 87.890625, 87.890625, 43.9453125, 
					 21.9726563, 21.9726563, 10.9863281, 5.49316406, 5.49316406, 2.74658203, 1.37329102, 1.37329102,
					 0.68664551, 0.34332275, 0.34332275, 0.17166138, 0.08583069, 0.02145767, 30000, 7500, 1875, 937.5,
					 468.75, 234.375, 117.1875, 117.1875, 58.59375, 29.296875, 29.296875, 14.6484375, 7.32421875,
					 7.32421875, 3.66210938, 1.83105469, 1.83105469, 0.91552734, 0.45776367, 0.45776367, 0.22888184,
					 0.11444092, 0.02861023)
}

## Calculate ERCC attomoles added per 1e5 cells
ERCC_Mix$attomoles_added = (ERCC_Mix$attomoles_per_ul)/(dilution)

## Calculate ERCC Moles
ERCC_Mix$moles = (ERCC_Mix$attomoles_added)/(1e18)

## Calculate ERCC molecules
ERCC_Mix$molecules =  (ERCC_Mix$moles)*(6.022e23)

## Calculate ERCC molecules per cell
ERCC_Mix$MPC = (ERCC_Mix$molecules)/(1e5)

## Calculate the log-base-2 of the MPC values
ERCC_Mix$log2MPC = log2(ERCC_Mix$MPC + 1)

## Merge ERCC log-base-2 MPC values and RPKM values for ERCC transcripts
StdCurve = merge(x = ERCC_Mix[,c(2,8)], y =  ERCC_abundance, by = "gene_id", all.y = TRUE)

## Calculate the log2(MPC) for each non-ERCC gene using x = (y-b)/m formula
## Experimental ERCC RPKM on y-axis and known ERCC MPC on x-axis
## Calculate the MPC using exponential 2x of log2(MPC) subtracted by 1

## Fit linear regression model
ERCC_lm = lm(log2RPKM ~ log2MPC, data = StdCurve)

## Get y-intercept
y_int = ERCC_lm$coefficients[[1]]

## Get slope
slope = ERCC_lm$coefficients[[2]]

## Calculate MPC for non-ERCC genes
abundance$MPC = round(((2^((abundance$log2RPKM-y_int)/slope))-1),2)

## Remove log2RPKM column
abundance = abundance[,-6]

## Calculate estimated counts with ERCC normalized results
abundance$expected_counts = round(((abundance$MPC)*(abundance$length/10^3))*(mapped_reads/10^6))

## Keep gene_id, expected_counts, and MPC
abundance = abundance[,c(-2:-3, -5)]

## Export ERCC normalized values as 
write.table(abundance, file=output_path, row.names = FALSE, quote = FALSE, sep = "\t")


## Differential Gene Expression Analsis with and/or without ERCC Spike-In RNA normalization

## Import libraries (DESeq2, tidyverse)
cat("Loading DESeq2 and tidyverse packages\n")
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))

cat("\nBeginning conventional Differential Gene Expression\n")
## Get inputs <UTS output directory> <uniquely mapped reads> <ERCC enabled/disabled> <ERCC Mix>
args = commandArgs(trailingOnly = TRUE)

## Set variables
counts_dir = args[1]
base_dir = dirname(counts_dir)
mapped_reads = as.numeric(unlist(strsplit(args[2], ' ')))
ERCC = args[3]
DGE_MIX = args[4]
sampleinfo = read.delim(args[5])
mapped_reads_names = unlist(strsplit(args[6], ' '))

## Create output directory
dir.create(file.path(base_dir, "DGE_output"), showWarnings = FALSE)

## Define location of output directory
DGE_output = file.path(base_dir, "DGE_output")

## Create data frame containing all samples and respective factors
study = sampleinfo

## Create data frame containing all samples and respective factors
study = sampleinfo[2]
rownames(study)=sampleinfo[[1]]

suppressWarnings( if (dim(study) >= 2){
	  group = apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
	  group = study[,1]
} )
group_names = paste0("(",group,")",sep = "")
group = make.names(group)
names(group) = group_names
rm(group_names)

## Format contrasts table
contrasts = combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names = combn(levels(factor(names(group))),2)
# format combinations for output table files names
contrast.names = c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),
		                      paste(contrast.names[2,],contrast.names[1,],sep = "v"))
contrasts = cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) = contrast.names
rm(contrast.names) 

## Import UTS output files as a list
files = grep(list.files(file.path(counts_dir), pattern = ".genes.results", full.names = TRUE), 
	      pattern = "ERCCnorm_", invert = TRUE, value = TRUE)

## Reorder files list according to sampleinfo
temp_list = as.list(sampleinfo[[1]])
temp_files = c()
for (i in as.character(temp_list)) {
	  temp_files = append(temp_files, c(files[grep(i, files)]))
}
files = temp_files
rm(temp_list, temp_files)

## Name Samples
names(files) = paste0(sampleinfo[[2]], "_", sampleinfo[[3]])

## Name uniquely mapped reads
names(mapped_reads) = paste0(mapped_reads_names)

## Reorder uniquely mapped reads according to sampleinfo
temp_list = as.list(sampleinfo[[1]])
temp_mapped_reads = c()
for (i in as.character(temp_list)) {
	temp_mapped_reads = append(temp_mapped_reads, c(mapped_reads[grep(i, names(mapped_reads))]))
}
mapped_reads = temp_mapped_reads
rm(temp_list, temp_mapped_reads)

## Import files as a matrix
df_list = suppressMessages(files %>%
	lapply(read_table))

for (i in 2:length(df_list)) {
	if (i == 2) {
		raw_counts = merge(as.data.frame(df_list[1])[,c(1,4)], as.data.frame(df_list[2])[,c(1,4)], all = TRUE, by = 1)
	} else {
		raw_counts = merge(raw_counts, as.data.frame(df_list[i])[,c(1,4)], all = TRUE, by = 1)
	}
}
rownames(raw_counts) = raw_counts[[1]]
raw_counts = raw_counts[,-1]
colnames(raw_counts) = c(names(df_list))
raw_counts[is.na(raw_counts)] = 0

## Make DESeqDataSet object
## Create data frame defining which group each sample belongs to
sampleTable = data.frame(condition=factor(group))
rownames(sampleTable) = colnames(raw_counts)

dds = DESeqDataSetFromMatrix(round(raw_counts), sampleTable, ~condition)

## Filter out genes with counts of less than 10 in all samples
keep = rowSums(counts(dds)) > 10
dds = dds[keep,]

## Generate a DESeqDataSet object without ERCC genes
dds_1 = dds

## Perform DESeq analysis without considering ERCC genes
size_factors = mapped_reads/mean(mapped_reads)
sizeFactors(dds_1) = size_factors
dds_1 = estimateDispersions(dds_1)
dds_1 = nbinomWaldTest(dds_1)
cat("\n")

## Get normalized counts
normCounts = as.data.frame(counts(dds_1, normalized=TRUE))

## Add 1 to all counts to avoid issues with log transformation
normCounts = normCounts + 1

## Generate F statistic p-value (similar to ANOVA p-value) using 
## DESeq2 likelihood ratio test (LRT) design
dds_1_lrt = DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt = results(dds_1_lrt)
cat("\n")

## Create reduced output table
reduced_output_table_1 = normCounts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts)[2]){
	  res_1 = results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
  res_1 = as.data.frame(res_1@listData)[,c(2,5,6)]
    colnames(res_1) =c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
    reduced_output_table_1 = cbind(reduced_output_table_1,res_1)
      rm(res_1)
}

## Generate and add all sample mean column to the normalized counts table
reduced_output_table_1$All.mean = rowMeans(normCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the (non-ERCC) normalized counts table
reduced_output_table_1$All.stdev = rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the (non-ERCC) normalized counts table
reduced_output_table_1$LRT.p.value = res_1_lrt@listData$padj

## Generate and add group mean and stdev columns to the (non-ERCC) normalized counts table
tcounts = as.data.frame(t(normCounts))
tcounts$group = group
group_means = as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means = group_means[-c(1),]
colnames(group_means) = paste0("Group.Mean_",levels(factor(names(group))))
group_stdev = as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev = group_stdev[-c(1),]
colnames(group_stdev) = paste0("Group.Stdev_",levels(factor(names(group))))

reduced_output_table_1 = cbind(reduced_output_table_1,group_means)
reduced_output_table_1 = cbind(reduced_output_table_1,group_stdev)
rm(group_stdev,group_means,tcounts)

## Export non-ERCC normalized DGE tables
write.table(reduced_output_table_1,file.path(DGE_output, "differential_expression"), row.names = TRUE, quote = FALSE, sep = "\t")

####################################### ERCC Normalized DGE ##########################################

if (ERCC == "ENABLED") {
  
  cat("\nBeginning ERCC Spike-In RNA normalized Differential Gene Expression\n")

  ## dds_2 will be used to generate data with considering ERCC genes
  dds_2 = dds

  ## Prepare filtered data to be normalized with ERCC genes
  ## Create list of rows containing ERCC genes to use for ERCC-normalization
  ## Note: ERCC genes should be the same concentration in all samples
  
  ## Create list of rows containing ERCC group B genes to use for ERCC-normalization
  ## Else, use all genes. Using group B is useful if Mix 1 and Mix 2 were used.
  ## Group B genes have same concentrations in both Mix 1 and Mix 2
  if (DGE_MIX == "group_b") {
	  ercc_rows = grep("ERCC-00096|ERCC-00171|ERCC-00009|ERCC-00042|ERCC-00060|ERCC-00035|ERCC-00025|ERCC-00051|ERCC-00053|ERCC-00148|ERCC-00126|ERCC-00034|ERCC-00150|ERCC-00067|ERCC-00031|ERCC-00109|ERCC-00073|ERCC-00158|ERCC-00104|ERCC-00142|ERCC-00138|ERCC-00117|ERCC-00075",rownames(dds_2))
    } else {
	    ercc_rows = grep("ERCC-", rownames(dds_2))
  }

    ## Identify and list sampels that do not contain counts for ERCC genes
    exit = function() { invokeRestart("abort") } # set exit function
    ercc_dds =  dds[ercc_rows,]  
    check_vector = (colSums(counts(ercc_dds)) == 0)
    check_vector = match("TRUE", check_vector)
    check_vector[is.na(check_vector)] = 0
    
    if (check_vector == 1) {
	    cat("Some samples do not have detectable ERCC spike-ins: ", colnames(ercc_dds[,colSums(counts(ercc_dds)) == 0]), sep="\n")
    }
    
    ## Perform DESeq analysis with considering ERC genes
    dds_2 = estimateSizeFactors(dds_2, controlGenes=ercc_rows)
    dds_2 = dds_2[-c(ercc_rows),] # remove ERCCs from counts table after normalization
    dds_2 = estimateDispersions(dds_2)
    dds_2 = nbinomWaldTest(dds_2)
    cat("\n")
    
    ## Get normalized counts
    ERCCnormCounts = as.data.frame(counts(dds_2, normalized=TRUE))
    
    ## Add 1 to all counts to avoid issues with log  transformation
    ERCCnormCounts = ERCCnormCounts + 1
    
    ## Generate F statistic p-value (similar to ANOVA p-value) using 
    ## DESeq2 likelihood ratio test (LRT) design
    dds_2_lrt = DESeq(dds_2, test = "LRT", reduced = ~ 1)
    res_2_lrt = results(dds_2_lrt)
    
    ## Create reduced output table
    reduced_output_table_2 = ERCCnormCounts
	      
    ## Iterate through Wald Tests to generate pairwise comparisons of all groups
    for (i in 1:dim(contrasts)[2]){
	    res_2 = results(dds_2, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	    res_2 = as.data.frame(res_2@listData)[,c(2,5,6)]
	    colnames(res_2)=c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
	    reduced_output_table_2 = cbind(reduced_output_table_2,res_2)
	    rm(res_2)
    }
    
    ## Generate and add all sample mean column to the (ERCC) normalized counts table
    reduced_output_table_2$All.mean = rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)
    
    ## Generate and add all sample stdev column to the (ERCC) normalized counts table
    reduced_output_table_2$All.stdev = rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)
    
    ## Add F statistic p-value (similar to ANOVA p-value) column to the (ERCC) normalized counts table
    reduced_output_table_2$LRT.p.value = res_2_lrt@listData$padj
    
    ## Generate and add group mean and stdev columns to the ERCC-normalized counts table
    tcounts = as.data.frame(t(ERCCnormCounts))
    tcounts$group = group
    group_means = as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
    group_means = group_means[-c(1),]
    colnames(group_means) = paste0("Group.Mean_",levels(factor(names(group))))
    group_stdev = as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
    group_stdev = group_stdev[-c(1),]
    colnames(group_stdev) = paste0("Group.Stdev_",levels(factor(names(group))))
			  
    reduced_output_table_2 = cbind(reduced_output_table_2,group_means)
    reduced_output_table_2 = cbind(reduced_output_table_2,group_stdev)
    rm(group_stdev,group_means,tcounts)
    
    ## Export ERCC normalized DGE tables
    write.table(reduced_output_table_2,file.path(DGE_output, "ERCCnorm_differential_expression"), row.names = TRUE, quote = FALSE, sep = "\t")    
}


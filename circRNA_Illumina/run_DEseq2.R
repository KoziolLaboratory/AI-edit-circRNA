library(tidyverse)
library(DESeq2)

dir <- "/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciriquant/bsj_files"

KO_filenames <- list.files(dir, pattern = "Adar2KO", full.names = TRUE)
WT_filenames <- list.files(dir, pattern = "WT", full.names = TRUE)

filenames <- c(KO_filenames, WT_filenames)
counts <- map(filenames, function(x) read.table(x, sep = "\t", stringsAsFactors = FALSE, header = TRUE))
names(counts) <- basename(filenames)
bsj_count_df <- map_dfc(counts, "bsj")
rownames(count_df) <- counts[[1]]$Geneid
colnames(count_df) <- gsub("\\..*", "", colnames(count_df))

# Now the coldata
# It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
coldata <- data.frame(condition = c(rep("Adar2KO", 3), rep("WT", 3)))
rownames(coldata) <- colnames(count_df)
all(rownames(coldata) == colnames(count_df))

dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = coldata,
                              design = ~ condition)

# Pre-filtering
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted.
# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# By default, R will choose a reference level for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels. To solve this, you can explicitly set the factors levels.

dds$condition <- factor(dds$condition, levels = c("Adar2KO","WT"))

dds <- DESeq(dds)
res_dds <- results(dds, contrast = c("condition", "WT", "Adar2KO"))
summary(res_dds)

# How many genes with adjusted p-values < 0.1 are there?
resLFC <- lfcShrink(dds, coef = "condition_WT_vs_Adar2KO", type = "apeglm")

# We can order our results table by the smallest p value:
resOrdered <- res_dds[order(res_dds$pvalue),]
sum(res_dds$padj < 0.1, na.rm=TRUE)

write.csv(resOrdered, file = file.path(dir, "Adar2KO_DEA.csv"))

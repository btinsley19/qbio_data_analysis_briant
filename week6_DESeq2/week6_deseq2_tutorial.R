#week6 DESeq2

library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)

BiocManager::install("DESeq2")
library(DESeq2)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method

# only need to download the data once! Comment this out once you have completed it once
#GDCdownload(query) 

sum_exp <- GDCprepare(query)

#make a copy of clinical data
patient_data <- colData(sum_exp)
#make a copy of counts data
counts <- assays(sum_exp)$"HTSeq - Counts"

#convert to dataframe
counts = data.frame(counts)

sum(is.na(patient_data$age_at_diagnosis))
#4 patients have NA for age

#remove NA data using boolean indexing
bool_age_na = is.na(patient_data$age_at_diagnosis)
patient_data <- patient_data[!bool_age_na,]
counts <- counts[,!bool_age_na]

patient_data$age_at_diagnosis = patient_data$age_at_diagnosis / 365
patient_data$age_category = ifelse(patient_data$age_at_diagnosis >= 50, "Old", "Young")

#convert age category to a factor
patient_data$age_category = factor(patient_data$age_category, levels = c("Young", "Old"))

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
#rownames(counts)
#rowData(sum_exp)$external_gene_name

counts_row_sums = rowSums(counts)
low_counts_mask = (counts_row_sums >= 10)
sum(low_counts_mask == FALSE)
#4817 genes have low counts

#remove low expressed genes from counts
counts <- counts[low_counts_mask,]


#time for analysis
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "Young", "Old"))

head(results)

row_order = order(results$padj)
results = results[row_order,]
head(results)
#ENSG00000202198 is significantly more expressed in older patients
#7SK RNA is a promoter / enhancer for rna expression

log2FoldChange_threshold = -1
padj_threshold = 0.05

genes_log2 <- results[results$log2FoldChange < log2FoldChange_threshold | results$log2FoldChange > abs(log2FoldChange_threshold),]
#genes_padj <- results[results$padj < padj_threshold,]



#this doesnt work
important_genes <- results[(results$log2FoldChange < log2FoldChange_threshold | 
                             results$log2FoldChange > abs(log2FoldChange_threshold)) & results$padj < padj_threshold,]



library(ggplot2)

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")

volcano_plot

write.csv(x = results,
          file = "/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week6_DESeq2/results.csv",
          row.names = FALSE)

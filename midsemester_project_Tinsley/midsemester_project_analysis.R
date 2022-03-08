#Brian Tinsley
#Midsemester Project

clinic <- data.table::fread("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/analysis_data/clinic.csv",
                            data.table = F)

### CLINICAL ANALYSIS ###

#look at columns of clinical data and choose one that is interesting
colnames(clinic)

#interested in that lymphatic columns
clinic$lymph_node_examined_count
clinic$lymphatic_invasion
clinic$primary_lymph_node_presentation_assessment #mostly yes, seems irrelevant for this study
clinic$number_of_lymphnodes_positive_by_he #curious what he and ihc are... google it
#he = 
clinic$number_of_lymphnodes_positive_by_ihc #this one has many NA's, lets ignore it
#ihc =Immunohistochemical

#Question: Are genes expressed differently in patients with/without lymphatic invasion?
#          Does the number of lymph nodes examined effect survival rate?

dim(clinic)
#there are 524 patients

#creating age category column for further analysis
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis > 50, "Old", "Young")

#clean up NA's and blank data
#new_clinic <- clinic[!is.na(clinic$lymph_node_examined_count),]
final_clinic <-new_clinic[new_clinic$lymphatic_invasion != "",]
#removed 65 patients from the data

sum(final_clinic$lymphatic_invasion == "YES") #188
sum(final_clinic$lymphatic_invasion == "NO") #271

youngyes <- sum(final_clinic$lymphatic_invasion =="YES" & final_clinic$age_category=="Young")
oldyes <- sum(final_clinic$lymphatic_invasion == "YES" & final_clinic$age_category=="Old")

data <- c(youngyes, oldyes)
#do older/younger patients experience lymphatic invasion more?
barplot(data, names=c("Young", "Old"), main="Patients with Lymphatic Invasion by Age")


#if there are more lymph nodes examined, do they have 
mean(final_clinic$lymph_node_examined_count)
median(final_clinic$lymph_node_examined_count)


#survival analysis time
library(survival)
library(survminer)

final_clinic$days_to_death <- ifelse(is.na(final_clinic$days_to_death), final_clinic$days_to_last_follow_up, final_clinic$days_to_death)
final_clinic$death_event <- ifelse(final_clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = final_clinic$days_to_death, 
                    event = final_clinic$death_event)

# We then create a fit object
invasion_fit <- surv_fit( surv_object ~ final_clinic$lymphatic_invasion, data = final_clinic )

survplot = ggsurvplot(invasion_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p




### RNA SEQ ANALYSIS ###
library(TCGAbiolinks)
library(SummarizedExperiment)
#load in the rna seq data
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
#only do this once
#GDCdownload(query)
sum_exp <- GDCprepare(query)
counts = assays(sum_exp)$"HTSeq - Counts"

colnames(colData(sum_exp))
colData(sum_exp)$age_at_index
colData(sum_exp)$paper_lymphatic_invasion_present
colData(sum_exp)$paper_lymphnodes_examined
colData(sum_exp)$paper_number_of_lymphnodes_examined
colData(sum_exp)$paper_lymphnode_pathologic_spread
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_index >= 50, "Old", "Young")
colData(sum_exp)$age_category

#choosing two genes to do analysis with, TP53 and APC
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "TP53")
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask] 
ensembl_geneA

geneB_id_mask = (rowData(sum_exp)$external_gene_name == "APC")
sum(geneB_id_mask) 
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask] 
ensembl_geneB

min(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ])
max(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ])

summary(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ])
summary(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneB, ])

plot(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneB, ],
     xlab = "TP53", ylab = "APC")

colData(sum_exp)$paper_lymphatic_invasion_present

bool_lymph_na = (is.na(colData(sum_exp)$paper_lymphatic_invasion_present) |
                   (colData(sum_exp)$paper_lymphatic_invasion_present == "NA") )
lymph_no_NAs <- colData(sum_exp)$paper_lymphatic_invasion_present[!bool_lymph_na]
length(lymph_no_NAs) #191 patients

gene_countsA = assays(sum_exp)$"HTSeq - Counts"[geneA_id_mask, !bool_lymph_na]
length(gene_countsA)
boxplot(gene_countsA ~ lymph_no_NAs, xlab = "Lymphatic Invasion (Yes/No)", ylab = "Counts for TP53 Gene", 
        main = "Counts for TP53 Gene for Patients with/without Lymphatic Invasion")

gene_countsB = assays(sum_exp)$"HTSeq - Counts"[geneB_id_mask, !bool_lymph_na]
length(gene_countsB)
boxplot(gene_countsB ~ lymph_no_NAs, xlab = "Lymphatic Invasion (Yes/No)", ylab = "Counts for APC Gene",
        main = "Counts for APC Gene for Patients with/without Lymphatic Invasion")





### DESEQ ANALYSIS ###
library(DESeq2)

#make a copy of clinical data
patient_data <- colData(sum_exp)
#make a copy of counts data
counts <- assays(sum_exp)$"HTSeq - Counts"

#convert to dataframe
counts = data.frame(counts)

bool_l <- (is.na(colData(sum_exp)$paper_lymphatic_invasion_present) |
             (colData(sum_exp)$paper_lymphatic_invasion_present == "NA") )
patient_data <- patient_data[!bool_l,]
patient_data$paper_lymphatic_invasion_present = factor(patient_data$paper_lymphatic_invasion_present,
                                                       levels=c("No", "Yes"))
counts <- counts[,!bool_l]

counts_row_sums = rowSums(counts)
low_counts_mask = (counts_row_sums >= 10)
sum(low_counts_mask == FALSE)
#8752 genes have low counts

#remove low expressed genes from counts
counts <- counts[low_counts_mask,]

#time to do analysis
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~paper_lymphatic_invasion_present)

dds_obj = DESeq(dds) #this takes awhile
resultsNames(dds_obj)  # see what comparisons got run
results = results(dds_obj, format = "DataFrame", contrast = c("paper_lymphatic_invasion_present", "NO", "YES"))

head(results)
row_order = order(results$padj)
results = results[row_order,]
head(results)

log2FoldChange_threshold = -1
padj_threshold = 0.05

genes_log2 <- results[results$log2FoldChange < log2FoldChange_threshold | results$log2FoldChange > abs(log2FoldChange_threshold),]

library(ggplot2)

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in NO lymph. invasion",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in NO lymph. invasion", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (NO/YES lymphatic invasion)",
       y = "-log10 Adjusted p-value")

volcano_plot




### MAF ANALYSIS ### 

library(maftools)

colnames(final_clinic)[ colnames(final_clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = final_clinic, 
                       isTCGA = TRUE)

maf_object
str(maf_object)

maf_object@data
maf_object@clinical.data$lymphatic_invasion
#both data sets have tumor sample barcode

oncoplot(maf = maf_object,
         top = 10)

clinic_maf <- maf_object@clinical.data

yes_patient_ids <- clinic_maf$Tumor_Sample_Barcode[clinic_maf$lymphatic_invasion=="YES"]
no_patient_ids <- clinic_maf$Tumor_Sample_Barcode[clinic_maf$lymphatic_invasion=="NO"]
yes_maf = subsetMaf(maf = maf_object,
                      tsb = yes_patient_ids)
no_maf = subsetMaf(maf = maf_object,
                    tsb = no_patient_ids)
coOncoplot(m1 = yes_maf, 
           m2 = no_maf, 
           m1Name = "YES lymphatic invasion", 
           m2Name = "NO lymphatic invasion")

lollipopPlot(maf=maf_object, gene="TP53")
lollipopPlot(maf=maf_object, gene="APC")

lollipopPlot2(m1 = yes_maf, 
              m2 = no_maf, 
              m1_name = "YES",
              m2_name = "NO",
              gene = "TP53")

lollipopPlot2(m1 = yes_maf, 
              m2 = no_maf, 
              m1_name = "YES",
              m2_name = "NO",
              gene = "APC")

# geneA_maf <- subsetMaf(maf = maf_object,
#                        genes = "TP53")
# 
# geneB_maf <- subsetMaf(maf = maf_object, 
#                        genes = "APC")
# 
# length(geneA_maf@data$Tumor_Sample_Barcode)
# length(geneB_maf@data$Tumor_Sample_Barcode)
# 
# mut_bc_geneA = geneA_maf@data$Tumor_Sample_Barcode
# mut_bc_geneB = geneB_maf@data$Tumor_Sample_Barcode
# 
# num_mut_geneA = length(mut_bc_geneA)
# num_mut_geneB = length(mut_bc_geneB)
# 
# mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
# num_mut_geneAB = length(mut_bc_geneAB)
# 
# num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
# num_mut_geneB_only = num_mut_geneB - num_mut_geneAB
# 
# num_mut_neither = 397 - (num_mut_geneA_only + num_mut_geneB_only + num_mut_geneAB)

#BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)
library(ggplot2)

clinic <- data.table::fread("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/analysis_data/clinic.csv",
                            data.table = F)

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

length(colnames(clinic))
length(colnames(clinic) == "bcr_patient_barcode")
colnames(clinic) == "bcr_patient_barcode"
typeof(colnames(clinic) == "bcr_patient_barcode")


mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#if we want to read in the maf data again we can find the csv
#maf_dataframe = data.table::fread("/Users/btinsley/Desktop/qbio_student_group/
#                                  qbio_data_analysis_briant/analysis_data/TCGA.COAD.
#                                 mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
#                                  data.table = F)
#maf_object2 <- read.maf(maf= maf_dataframe, clinicalData = clinic, isTCGA = TRUE)

maf_object
str(maf_object)

maf_object@data
maf_object@clinical.data
#both data sets have tumor sample barcode

oncoplot(maf = maf_object,
         top = 10)

ggsave("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week7_MAF/oncoplot.png")

#the KRAS gene is the 4th highest mutated. The KRAS gene provides instructions 
#for making a protein called K-Ras that is part of a signaling pathway 
#known as the RAS/MAPK pathway. The protein relays signals from outside the cell to the cell's nucleus

clinic <- maf_object@clinical.data

young_patient_ids <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$age_at_initial_pathologic_diagnosis<50]
old_patient_ids <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$age_at_initial_pathologic_diagnosis>=50]
young_maf = subsetMaf(maf = maf_object,
                      tsb = young_patient_ids)
old_maf = subsetMaf(maf = maf_object,
                      tsb = old_patient_ids)
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young", 
           m2Name = "Old")
ggsave("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week7_MAF/young_old_oncoplot.png")


lollipopPlot(maf=maf_object, gene="KRAS")
ggsave("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week7_MAF/lollipop.png")

# pick a gene to look at!
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "Young",
              m2_name = "Old",
              gene = "KRAS")
ggsave("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week7_MAF/young_old_lollipop.png")

#seems to be mutated about the same in young and old

#b=7, c=2, d=35, e=37 , f=42
#if you have no mutation in gene A, you are very lkely to have no mutation in gene B

# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

geneA_maf
geneB_maf

length(geneA_maf@data$Tumor_Sample_Barcode)
length(geneB_maf@data$Tumor_Sample_Barcode)

mut_bc_geneA = geneA_maf@data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@data$Tumor_Sample_Barcode

num_mut_geneA = length(mut_bc_geneA)
num_mut_geneB = length(mut_bc_geneB)

mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)

num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB

num_mut_neither = 397 - (num_mut_geneA_only + num_mut_geneB_only + num_mut_geneAB)

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_mut_neither), 
                       nrow=2)

contig_table

fe_results <- fisher.test(contig_table)
fe_results
#these results yield a low p value so there is a correlation betweent the two mutations



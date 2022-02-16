library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method

# only need to download the data once! Comment this out once you have completed it once
#GDCdownload(query) 

sum_exp <- GDCprepare(query)

counts = assays(sum_exp)$"HTSeq - Counts"
counts
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]
head(rowData(sum_exp))
head(colData(sum_exp))

colData(sum_exp)[1:5, 25:29]

dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")
#the rows of coldata are same size as the columns of heseq-counts
#the rows of row data are the same as rows of htseq-counts

str(colData(sum_exp))
head(colData(sum_exp))
#rows of coldata are columns of htseq counts, the columns are data about each patient
#str(rowData(sum_exp))
head(rowData(sum_exp))

colnames(colData(sum_exp))
(colData(sum_exp)$age_at_diagnosis)
colData(sum_exp)$age_at_diagnosis[1:10]
#units are days

colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis/365
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis >= 50, "Old", "Young")
colData(sum_exp)$age_category

head(rowData(sum_exp))
dim(rowData(sum_exp))
#each row is one gene and it contains the different names for that gene

"TP53" %in% rowData(sum_exp)$external_gene_name
"APC" %in% rowData(sum_exp)$external_gene_name

assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]


geneA_id_mask = (rowData(sum_exp)$external_gene_name == "TP53")
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask] #fill in the dotted lines. 
ensembl_geneA

geneB_id_mask = (rowData(sum_exp)$external_gene_name == "APC")
sum(geneB_id_mask) 
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask] #fill in the dotted lines. 
ensembl_geneB

min(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ])
max(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ])

summary(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneB, ])

plot(assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ ensembl_geneB, ],
     xlab = "TP53", ylab = "APC")

sum(is.na(colData(sum_exp)$age_category))
#colData(sum_exp)$age_category

bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na

age_cat_no_NAs <- colData(sum_exp)$age_category[!bool_age_na]
length(age_cat_no_NAs)

dim( colData(sum_exp) )[1] == length(age_cat_no_NAs)+num_na
#columns of colData is the same as the columns of new dataframe

dim(assays(sum_exp)$"HTSeq - Counts")
#there are 521 patients

dim(age_cat_no_NAs)
#there are 4 less patients here because we removed the NA's

identical( rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")  )
#returns true

gene_counts = assays(sum_exp)$"HTSeq - Counts"[geneA_id_mask, !bool_age_na]
length(gene_counts)

boxplot(gene_counts ~ age_cat_no_NAs, xlab = "Age Category (Young/Old)", ylab = "Counts for TP53 Gene")


#1. How would you access the HTSeq - Counts data frame? What data do the rows represent? What about the columns?
assays(sum_exp)$"HTSeq - Counts"
#rows represent genes
#columns represent the counts of that gene for different patients 

#2. How would you access the data frame that stores more information about the rows in the counts data frame? 
#Specifically, what information is being stored? How does this data frame relate to the rows of the counts data frame?
rowData(sum_exp)
#row data has different names for the gene in each row. The rows in this data frame correspond to the rows in counts

#3.How would you access the data frame that stores more information about the columns in the counts data frame? 
#Specifically, what information is being stored? How does this data frame relate to the columns of the counts data frame?
colData(sum_exp)
#data about the patients is being stored. The rows of coldata are the columns in counts data frame
#clinical data is stored here, so there are a bunch of columns for each patient with different info


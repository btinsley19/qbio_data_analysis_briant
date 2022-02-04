if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
if(!require(TCGAbiolinks)) {BiocManager::install("TCGAbiolinks")}
library(TCGAbiolinks)
# You will use variations of this command for every data type
clin_query <- GDCquery(project = "TCGA-COAD",
                       data.category = "Clinical",
                       file.type = "xml")
# ONLY RUN THIS ONCE
#GDCdownload(clin_query)

# Load the query into R
clinic <- GDCprepare_clinic(clin_query,clinical.info = "patient")
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

plot(clinic$weight, clinic$age, xlab="weight", ylab="age")
unique(clinic$race_list)
par(mar=c(10,1,1,1))
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list, las=2)

sum(clinic$race_list == "")
clinic$race_list <- ifelse(clinic$race_list == "", "No Data", as.character(clinic$race_list))

min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
summary(clinic$age_at_initial_pathologic_diagnosis)

dim(clinic)
sum(clinic$age_at_initial_pathologic_diagnosis < 50)
sum(clinic$age_at_initial_pathologic_diagnosis >= 50)

young_patient_ids <- clinic$bcr_patient_barcode[clinic$age_at_initial_pathologic_diagnosis < 50]
old_patient_ids <- clinic$bcr_patient_barcode[clinic$age_at_initial_pathologic_diagnosis >= 50]

clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "Young", "Old")

clinic[1,1] #This is the top left entry of the dataframe. R has "one-based indexing"
clinic[1,] #this is the first row of the dataframe (1 patients information)
clinic[2:5,] #this is the 2nd to 4th rows of the dataframe 
clinic[,3] #this is the 3rd colmn with all rows

young_clinic <- clinic[clinic$bcr_patient_barcode %in% young_patient_ids,]
old_clinic <- clinic[clinic$bcr_patient_barcode %in% old_patient_ids,]


write.csv(clinic, "/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week4_clinical/coad_clinical_data.csv", row.names = F)


install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)

clinic$days_to_death <- ifelse(is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death)
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
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

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("/Users/btinsley/Desktop/qbio_student_group/qbio_data_analysis_briant/week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)



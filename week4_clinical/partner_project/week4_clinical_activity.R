# Brian Tinsley

clinic <- read.csv("clinic.csv")

# 1. Define the following: categorical variable, discrete variable, continuous variable. Provide examples of each
# categorical variables -- can be non-numeric and the data fits into certain groups, and there is no specific ordering most of
# the time such as color
# discrete variables -- countable numeric variables within two values
# continuous variable -- infinite number of variables between two values (like data and time)


#2. 
colnames(clinic)
clinic$lymphatic_invasion
sum(clinic$lymphatic_invasion == "")

#3
# Lymphatic invasion is measured by taking a biopsy of tissues
# surrounding the tumor and examining under a microscope
# This variable is categorical


#4
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2605104/
# this article explains that the effect of LVI on prognosis 
# widely unknown, and that a need for criteria in evaluation of 
# lymphovascular invasion is necessary for better prognosis

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5542293/
# article explains that LVI is significantly associated with 
# poor prognosis in overall survival

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4807319/
# this article considers the numerous factors that affect lymph 
# node yields, including the impact of the surgical and histopathological 
# techniques, as well as potential future strategies, including the use of 
# evolving technologies to increase efficiency of lymph node data collection

#5
colnames(clinic)
clinic$lymph_node_examined_count
# this variable determines the number of lymph nodes that were examined
# in the biopsies 

#6
# Higher levels of lymph nodes counted yields higher chances of lymphatic invasion
# Lymphatic invasion yields lower levels of survival in patients
# more lymph nodes studied yields higher survival and better prognosis in patients

#7
# below, i plotted lymph node invasion (Yes/No) and the number of lymph nodes
# examined. The plot did not show much. The survival plots, however, show that
# lower lymph node count and lymphatic invasion both lead to poorer survivability

#clean up data from NA's
new_clinic <- clinic[!is.na(clinic$lymph_node_examined_count),]
mean(new_clinic$lymph_node_examined_count)
# turn continous variable into categorical column
new_clinic$lymph_count = ifelse(new_clinic$lymph_node_examined_count > 22, "Lot", "Little")
#clean up data for other variable with blank data
final_clinic <-new_clinic[new_clinic$lymphatic_invasion != "",]
#plot 
boxplot(final_clinic$lymph_node_examined_count ~ final_clinic$lymphatic_invasion)

#survival analysis time
#install packages
install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer)

final_clinic$days_to_death <- ifelse(is.na(final_clinic$days_to_death), final_clinic$days_to_last_follow_up, final_clinic$days_to_death)
final_clinic$death_event <- ifelse(final_clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = final_clinic$days_to_death, 
                    event = final_clinic$death_event)

# We then create a fit object
invasion_fit <- surv_fit( surv_object ~ final_clinic$lymphatic_invasion, data = final_clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
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



# now going to do it for lymph count

# We then create a fit object
lymph_count_fit <- surv_fit( surv_object ~ final_clinic$lymph_count, data = final_clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(lymph_count_fit, 
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


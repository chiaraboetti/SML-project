library(stringi)
library(ggplot2)
library(dplyr)
library(knitr)
library(smotefamily)
`%notin%` <- Negate(`%in%`)
#####################################################

df1 = read.csv("dataset/CRISPR_gene_dependency.csv")
df2 = read.csv("dataset/sample_info.csv")

# checking if label of df1 are all contained in df2
prod(rownames(df1) %in% rownames(df2))


############################################################
# ??? # (de-)Capitalize columns in df2 and filter it
############################################################

unique(df2$primary_disease)
df2$primary_disease = stri_trans_totitle(df2$primary_disease)
unique(df2$primary_disease)

unique(df2$sample_collection_site)
df2$sample_collection_site = chartr(old = "_", new = " ", 
                                    df2$sample_collection_site)
df2$sample_collection_site = stri_trans_totitle(df2$sample_collection_site)
unique(df2$sample_collection_site)

unique(df2$lineage)
df2$lineage = chartr(old = "_", new = " ", df2$lineage)
df2$lineage = stri_trans_totitle(df2$lineage)
unique(df2$lineage)

trunc_df2 = df2 %>%
  filter(DepMap_ID %in% df1$DepMap_ID)

############################################################
# ??? # Looking for weird obs labels
############################################################
which(trunc_df2$primary_disease == "Unknown")
which(trunc_df2$primary_disease == "")
which(trunc_df2$primary_disease == "Non-Cancerous")
# --> c(986, 996)
which(trunc_df2$primary_disease == "Immortalized")
which(trunc_df2$primary_disease == "Engineered")
# --> c(715, 1024, 1025, 1026, 1027, 1028)

# Investigating more about "Non-Cancerous" and "Engineered"
# Since also lineage can be used as Cancer type classifications,
# we can use it as double check:

trunc_df2 %>%
  filter(primary_disease %in% c("Non-Cancerous","Engineered")) %>% 
  select(primary_disease, Subtype, sample_collection_site,
         lineage, lineage_subtype) %>%
  kable()

# All these obs comes from Engineered lineage and have no indication
# about the subtype disease and the subtype lineage.
# This makes sense as engineered cells have been synthetically modified.
# Thus, two ways to proceed:
# AUT delete 8 obs out to 1032,
# AUT keep them all, but putting them into the classes wtr the sample
# collection site.

# Conclusion: we remove the 2 Non-Cancerous cells and keep
#             Engineered cells

############################################################
# ??? # Checking for NAs
############################################################
nas = data.frame(which(is.na(df1), arr.ind = TRUE))
nas = nas[order(nas$row), ]
nas = cbind(df1$DepMap_ID[nas$row], nas)
colnames(nas) = c("DepMap_ID", "row", "col") 

ggplot(nas, aes(x = DepMap_ID)) +
  geom_bar(aes(fill = DepMap_ID), stat = "count") +
  ggtitle("How many NA per Cell") +
  labs(x = "Cell DepMap IDs", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

# Let us see how many missing values per cell
NaCount = nas %>%
  group_by(row) %>%
  summarise(count = n())
NaCount = cbind(df1$DepMap_ID[unique(nas$row)], NaCount)
colnames(NaCount) = c("DepMap_ID", "Row_index", "Count")
kable(NaCount)
# REM: weird obs are c(715, 986, 996, 1024, 1025, 1026, 1027, 1028)
# --> No correspondence with Na

trunc_df2 %>%
  filter(DepMap_ID %in% nas$DepMap_ID) %>%
  select(DepMap_ID, primary_disease, sample_collection_site, lineage) %>%
  kable()
# Nothing particular, obs seem unrelated

# Since there are only ten obs which have NA values, we decide 
# to remove them all:
df1 = df1 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)
trunc_df2 = trunc_df2 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)

data.frame(which(is.na(df1), arr.ind = TRUE))
data.frame(which(is.na(trunc_df2$primary_disease), arr.ind = TRUE))

# Weird obs:
which(trunc_df2$primary_disease %in% c("Non-Cancerous", "Engineered"))
# --> c(707, 976, 986, 1014, 1015, 1016, 1017, 1018)


############################################################
# ??? # Check if df1 is a data.frame of probabilities
############################################################

## check if entries are in [0, 1]: ok!
prod((df1 %>% select(-DepMap_ID)) >= 0 & (df1 %>% select(-DepMap_ID)) <= 1)

## check of row sum
row.sum = rowSums(df1 %>% select(-DepMap_ID))

# We observe that the sum of the rows does not sum to 1. 
# That is not a problem despite dealing with probabilities: 
# in fact, the effect of a given gene on a certain cancer 
# does not exclude the effect of another gene on it.


############################################################
# ??? # Comparison btw primary_disease and sample_collection_site
############################################################

cat("n. diseases =", length(unique(trunc_df2$primary_disease)), 
    "VS. n. sites =", length(unique(trunc_df2$sample_collection_site)))
# more sites then disease: let us investigate more and 
# see if there is an evident substructure

CancerSiteCount = trunc_df2 %>%
  select(DepMap_ID, primary_disease, sample_collection_site) %>%
  rename(disease = primary_disease,
         site = sample_collection_site)

CancerSiteCount %>%
  group_by(disease, site) %>% 
  summarise(count = n()) %>%
  kable()
# Peculiarities. probably because of metastasis:
# for instance, some Brain cancer cells have been collected
# from the abdomen, Lung cancer cells from a variety of 
# different places

ggplot(CancerSiteCount, aes(x = disease)) +
  geom_histogram(aes(fill = site), stat = "count") +
  ggtitle("Number of cancer per site") +
  labs(x = "Type of cancer", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

# Conclusion: we retain that it is not a good idea
# looking at the sample sites for deciding how to
# group cancer types together
# Another idea could be using an unsupervised cluster
# algorithm


############################################################
# ??? # Count obs wrt primary disease
############################################################

ggplot(trunc_df2, aes(x = primary_disease)) +
  geom_histogram(aes(fill = primary_disease), stat = "count") +
  ggtitle("Number of cancer types") +
  labs(x = "Type of cancer", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

CancerCount = trunc_df2 %>%
  group_by(primary_disease) %>%
  summarise(count = n())
kable(CancerCount)

CancerCount %>%
  filter(primary_disease %in% c("Bile Duct Cancer", "Kidney Cancer", "Gastric Cancer",
                                "Gallbladder Cancer", "Esophageal Cancer", 
                                "Colon/Colorectal Cancer", "Liver Cancer")) %>%
  select(count) %>%
  sum()

CancerCount %>%
  filter(primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma")) %>%
  select(count) %>%
  sum()

############################################################
# ??? # Obtain a dataset for binary classification 
############################################################

# We decide which are the labels related to various types of
# cancer according to [https://www.cancer.gov/types/by-body-location]
# In particular, different classification problems:
# - multi-classification
# - lung vs all: fewer number of obs but homogeneous
#   (126 obs)
# - blood vs all: fewer number of obs but homogeneous
#   (109 obs)


############################################################
# ??? # Gastrointestinal cancers VS. All

trunc_df2 = trunc_df2 %>%
  mutate(isGastro = case_when (
    primary_disease == "Bile Duct Cancer" ~ 1,
    primary_disease == "Kidney Cancer" ~ 1,
    primary_disease == "Gastric Cancer" ~ 1,
    primary_disease == "Gallbladder Cancer" ~ 1,
    primary_disease == "Esophageal Cancer" ~ 1,
    primary_disease == "Colon/Colorectal Cancer" ~ 1,
    primary_disease == "Liver Cancer"~ 1,
    TRUE ~ 0 
  ))

# appending this column to the CRISP dataset to have labeled data
label = trunc_df2$isGastro
gastro_dataset = cbind(df1, label)

# and splitting the data into training and test 
set.seed(8675309)
n.train = floor(.80*dim(df1)[1])
training = sample(1:dim(df1)[1], size = n.train, replace = FALSE)

gastro.training = gastro_dataset[training, ]  
gastro.test = gastro_dataset[-training, ]

# producing csv file of training and test
# write.csv(gastro.training,"C:/Users/user/Desktop/SML/Project/gastro_training.csv", row.names = F)
# write.csv(gastro.test,"C:/Users/user/Desktop/SML/Project/gastro_test.csv", row.names = F)


############################################################
# ??? # Lung cancer VS. All

trunc_df2 = trunc_df2 %>%
  mutate(isLung = as.numeric(primary_disease == "Lung Cancer"))

# appending this column to the CRISP dataset to have labeled data
label = trunc_df2$isLung
lung_dataset = cbind(df1, label)

# and splitting the data into training and test 
set.seed(8675309)
n.train = floor(.80*dim(df1)[1])
training = sample(1:dim(df1)[1], size = n.train, replace = FALSE)

lung.training = lung_dataset[training, ]  
lung.test = lung_dataset[-training, ]

# producing csv file of training and test
# write.csv(lung.training,"C:/Users/user/Desktop/SML/Project/lung_training.csv", row.names = F)
# write.csv(lung.test,"C:/Users/user/Desktop/SML/Project/lung_test.csv", row.names = F)

# balancing the minority class using the SMOTE function 
length(which(lung.training$label == 1))/nrow(lung.training)

genData = SMOTE(lung.training[,-1], lung.training$label, dup_size = 1)
lung.training_balanced = genData$data

length(which(lung.training_balanced$label == 1))/nrow(lung.training_balanced)
# 192 positive obs over 917 (about 21%)

# write.csv(lung.training_balanced,"C:/Users/user/Desktop/SML/Project/lung_training_balanced.csv", 
#           row.names = F)


############################################################
# ??? # Blood cancer VS. All

trunc_df2 = trunc_df2 %>%
  mutate(isBlood = case_when (
    primary_disease == "Leukemia" ~ 1,
    primary_disease == "Lymphoma" ~ 1,
    primary_disease == "Myeloma" ~ 1,
    TRUE ~ 0 
  ))

# appending this column to the CRISP dataset to have labeled data
label = trunc_df2$isBlood
blood_dataset = cbind(df1, label)

# and splitting the data into training and test 
set.seed(8675309)
n.train = floor(.80*dim(df1)[1])
training = sample(1:dim(df1)[1], size = n.train, replace = FALSE)

blood.training = blood_dataset[training, ]  
blood.test = blood_dataset[-training, ]

# producing csv file of training and test
# write.csv(blood.training,"~/Documents/sds/sml/SML-project/dataset/blood_training.csv", row.names = F)
# write.csv(blood.test,"~/Documents/sds/sml/SML-project/dataset/blood_test.csv", row.names = F)



############################################################
# ??? # Obtain dataset for multiclass classification
############################################################

trunc_df2_bis = trunc_df2 %>%
  filter(primary_disease != "Non-Cancerous") %>%
  mutate(CancerType = case_when (
    primary_disease == "Eye Cancer" ~ 0,
    (primary_disease == "Engineered") & (sample_collection_site == "Eye") ~ 0,
    
    primary_disease == "Bile Duct Cancer" ~ 1,
    primary_disease == "Pancreatic Cancer" ~ 1,
    primary_disease == "Gastric Cancer" ~ 1,
    primary_disease == "Gallbladder Cancer" ~ 1,
    primary_disease == "Esophageal Cancer" ~ 1,
    primary_disease == "Colon/Colorectal Cancer" ~ 1,
    primary_disease == "Liver Cancer"~ 1,
    
    primary_disease == "Ovarian Cancer" ~ 2,
    primary_disease == "Cervical Cancer" ~ 2,
    primary_disease == "Endometrial/Uterine Cancer" ~ 2,
    primary_disease == "Teratoma" ~ 2,
    
    primary_disease == "Bone Cancer" ~ 3,
    primary_disease == "Skin Cancer" ~ 3,
    primary_disease == "Fibroblast" ~ 3,
    primary_disease == "Sarcoma" ~ 3,
    primary_disease == "Liposarcoma" ~ 3,
    
    primary_disease == "Brain Cancer" ~ 4,
    primary_disease == "Neuroblastoma" ~ 4,
    primary_disease == "Rhabdoid" ~ 4,
    
    primary_disease == "Breast Cancer" ~ 5,
    
    primary_disease == "Head And Neck Cancer" ~ 6,
    primary_disease == "Thyroid Cancer" ~ 6,
    
    primary_disease == "Leukemia"~ 7,
    primary_disease == "Lymphoma"~ 7,
    primary_disease == "Myeloma"~ 7,
    
    primary_disease == "Bladder Cancer" ~ 8,
    primary_disease == "Kidney Cancer" ~ 8,
    primary_disease == "Prostate Cancer" ~ 8,
    (primary_disease == "Engineered") & (sample_collection_site == "Kidney") ~ 8,
    
    primary_disease == "Lung Cancer" ~ 9,
    
    TRUE ~ 10 
  ))

label = trunc_df2_bis$CancerType

df1_bis = df1 %>%
  filter(DepMap_ID %in% trunc_df2_bis$DepMap_ID)
multiclass_dataset = cbind(df1_bis, label)

set.seed(8675309)
n.train = floor(.80*dim(df1_bis)[1])
training = sample(1:dim(df1_bis)[1], size = n.train, replace = FALSE)

multiclass.training = multiclass_dataset[training, ]  
multiclass.test = multiclass_dataset[-training, ]

# write.csv(multiclass.training,"C:/Users/user/Desktop/SML/Project/multiclass_training.csv", 
#           row.names = F)
# write.csv(multiclass.test,"C:/Users/user/Desktop/SML/Project/multiclass_test.csv", 
#           row.names = F)



ggplot(trunc_df2_bis, aes(x = as.factor(CancerType))) +
  geom_bar(aes(fill = primary_disease), width = 0.55, 
           position = position_dodge(width = 0.75)) +
  ggtitle("Cancers by Body Location/System") +
  labs(x = "Location", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Eye", "Gastrointestinal", "Gynecologic", "Musculoskeletal",
                              "Neurologic", "Breast", "Head-Neck", "Hematologic",
                              "Genitourinary", "Lung")) +
  scale_fill_discrete(name = "Primary disease")

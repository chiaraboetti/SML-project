library(stringi)
library(ggplot2)
library(dplyr)
library(knitr)
library(smotefamily)
`%notin%` <- Negate(`%in%`)
#####################################################

df1 = read.csv("CRISPR_gene_dependency.csv")
df2 = read.csv("sample_info.csv")

# checking if label of df1 are all contained in df2
prod(rownames(df1) %in% rownames(df2))


############################################################
# ♫ # (de-)Capitalize columns in df2 and filter it
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
# ♫ # Looking for weird obs labels
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


############################################################
# ♫ # Checking for NAs
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

trunc_df2 %>%
  filter(DepMap_ID %in% nas$DepMap_ID) %>%
  select(DepMap_ID, primary_disease, sample_collection_site, lineage) %>%
  kable()
# Nothing particular, obs seem unrelated

# Since there are only ten obs which have NA values, we decide 
#to remove them all:
df1 = df1 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)
trunc_df2 = trunc_df2 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)

data.frame(which(is.na(df1), arr.ind = TRUE))
data.frame(which(is.na(trunc_df2$primary_disease), arr.ind = TRUE))

# NaCount[NaCount$Count > 1000, ]$Row_index # to be deleted
# 
# new_df1 = df1 %>%
#   filter(DepMap_ID %notin% NaCount[NaCount$Count > 1000, ]$DepMap_ID)
# 
# # And let us replace NA with 0
# new_df1[is.na(new_df1)] = 0


############################################################
# ♫ # Check if df1 is a data.frame of probabilities
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
# ♫ # Comparison btw primary_disease and sample_collection_site
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
# Stranger things: probably it is not a good idea looking at the
# sample sites for the clusters

ggplot(CancerSiteCount, aes(x = disease)) +
  geom_histogram(aes(fill = site), stat = "count") +
  ggtitle("Number of cancer per site") +
  labs(x = "Type of cancer", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))
# Conclusion: just stick to the primary_disease for clusters 


############################################################
# ♫ # Count obs wrt primary disease
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
  sum(select(count))


############################################################
# ♫ # Obtain a dataset for binary classification 
############################################################

# We decide to try two binary classification problems:
# - gastrointestinal vs all: this class has the highest number 
#                            of obs but heterogeneous
#   (205 obs)
# - lung vs all: fewer number of obs but homogeneous
#   (126 obs)

############################################################
# Gastrointestinal cancers VS. All the others
# We decide which are the labels related to Gastrointestinal
# cancers according to [https://www.cancer.gov/types/by-body-location]

trunc_df2 = trunc_df2 %>%
  mutate(isGastro = case_when (
    primary_disease ==  "Bile Duct Cancer" ~ 1,
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
# ♫ # Lung cancer VS. All

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
genData = SMOTE(lung.training[,-1], lung.training$label, dup_size = 4)
lung.training_balanced = genData$data
length(which(lung.training_balanced$label == 1))/nrow(lung.training_balanced)
# (about 40%)

# write.csv(lung.training_balanced,"~/Documents/sds/sml/SML-project/dataset/lung_training_balanced.csv", 
#           row.names = F)


############################################################
# ♫ # Obtain dataset for 3/4 classes classification
############################################################

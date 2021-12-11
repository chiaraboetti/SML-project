# DATA MANIPULATION

library(stringi)
library(ggplot2)
library(dplyr)
library(knitr)
`%notin%` <- Negate(`%in%`)
#####################################################

## load the datasets
df1 = read.csv("dataset/CRISPR_gene_dependency.csv")
df2 = read.csv("dataset/sample_info.csv")

# check if label of df1 are all contained in df2
prod(rownames(df1) %in% rownames(df2))

# preparation of trunc_df2
df2$primary_disease = stri_trans_totitle(df2$primary_disease)
unique(df2$primary_disease)

df2$sample_collection_site = chartr(old = "_", new = " ", 
                                    df2$sample_collection_site)
df2$sample_collection_site = stri_trans_totitle(df2$sample_collection_site)
unique(df2$sample_collection_site)

df2$lineage = chartr(old = "_", new = " ", df2$lineage)
df2$lineage = stri_trans_totitle(df2$lineage)
unique(df2$lineage)

trunc_df2 = df2 %>%
  filter(DepMap_ID %in% df1$DepMap_ID)


## look for weird obs labels
which(trunc_df2$primary_disease == "Unknown")
which(trunc_df2$primary_disease == "")
which(trunc_df2$primary_disease == "Non-Cancerous")
# --> c(986, 996)
which(trunc_df2$primary_disease == "Immortalized")
which(trunc_df2$primary_disease == "Engineered")
# --> c(715, 1024, 1025, 1026, 1027, 1028)

# investigate "Non-Cancerous" and "Engineered"
trunc_df2 %>%
  filter(primary_disease %in% c("Non-Cancerous","Engineered")) %>% 
  select(primary_disease, Subtype, sample_collection_site,
         lineage, lineage_subtype) %>%
  kable()

# All these obs comes from Engineered lineage and have no indication about the
# subtype disease and the subtype lineage. This makes sense as engineered cells
# have been synthetically modified


## look for NAs
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

# Nothing particular, obs belong to different cancer and there are only 10
#   --> remove them all
df1 = df1 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)
trunc_df2 = trunc_df2 %>%
  filter(DepMap_ID %notin% NaCount$DepMap_ID)

data.frame(which(is.na(df1), arr.ind = TRUE))
data.frame(which(is.na(trunc_df2$primary_disease), arr.ind = TRUE))

which(trunc_df2$primary_disease %in% c("Non-Cancerous", "Engineered"))
# weird obs are c(707, 976, 986, 1014, 1015, 1016, 1017, 1018)


## check if df1 is a data.frame of probabilities
prod((df1 %>% select(-DepMap_ID)) >= 0 & (df1 %>% select(-DepMap_ID)) <= 1)

row.sum = rowSums(df1 %>% select(-DepMap_ID))
# Sum of the rows does not sum to 1, but not a problem. 
# In fact, the effect of a given gene on a certain cancer does not exclude
# the effect of another gene on it.


## comparison btw primary_disease and sample_collection_site
cat("n. diseases =", length(unique(trunc_df2$primary_disease)), 
    "VS. n. sites =", length(unique(trunc_df2$sample_collection_site)))

CancerSiteCount = trunc_df2 %>%
  select(DepMap_ID, primary_disease, sample_collection_site) %>%
  rename(disease = primary_disease,
         site = sample_collection_site)

CancerSiteCount %>%
  group_by(disease, site) %>% 
  summarise(count = n()) %>%
  kable()
# Metastasis:
# for instance, some Brain cancer cells have been collected from the abdomen,
# Lung cancer cells from a variety of different places

ggplot(CancerSiteCount, aes(x = disease)) +
  geom_histogram(aes(fill = site), stat = "count") +
  ggtitle("Number of cancer per site") +
  labs(x = "Type of cancer", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

# Conclusion: we retain that it is not a good idea looking at the sample sites
# for deciding how to group cancer types together 
# Another idea could be using an unsupervised cluster algorithm


## primary disease
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
## dataset for binary classifications 
# RMK: we decide which are the labels related to various types of cancer
# according to [https://www.cancer.gov/types/by-body-location]

## Gastrointestinal cancers VS. All
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
# write.csv(gastro.training, "dataset/gastro_training.csv", row.names = F)
# write.csv(gastro.test, "dataset/gastro_test.csv", row.names = F)


## Lung cancer VS. All
trunc_df2 = trunc_df2 %>%
  mutate(isLung = as.numeric(primary_disease == "Lung Cancer"))

label = trunc_df2$isLung
lung_dataset = cbind(df1, label)

# training and test sets 
set.seed(8675309)
n.train = floor(.80*dim(df1)[1])
training = sample(1:dim(df1)[1], size = n.train, replace = FALSE)

lung.training = lung_dataset[training, ]  
lung.test = lung_dataset[-training, ]

# write.csv(lung.training, "dataset/lung_training.csv", row.names = F)
# write.csv(lung.test, "dataset/lung_test.csv", row.names = F)


## Blood cancer VS. All
trunc_df2 = trunc_df2 %>%
  mutate(isBlood = case_when (
    primary_disease == "Leukemia" ~ 1,
    primary_disease == "Lymphoma" ~ 1,
    primary_disease == "Myeloma" ~ 1,
    TRUE ~ 0 
  ))

label = trunc_df2$isBlood
blood_dataset = cbind(df1, label)

# training and test sets
set.seed(8675309)
n.train = floor(.80*dim(df1)[1])
training = sample(1:dim(df1)[1], size = n.train, replace = FALSE)

blood.training = blood_dataset[training, ]  
blood.test = blood_dataset[-training, ]

# write.csv(blood.training, "dataset/blood_training.csv", row.names = F)
# write.csv(blood.test, "dataset/blood_test.csv", row.names = F)



############################################################
## dataset for multiclass classification

# Multiclass with NO weird obs
trunc_df2_bis = trunc_df2 %>%
  filter(primary_disease %notin% c("Non-Cancerous","Engineered", "Eye Cancer")) %>%
  mutate(CancerType = case_when (
    primary_disease == "Bile Duct Cancer" ~ 0,
    primary_disease == "Pancreatic Cancer" ~ 0,
    primary_disease == "Gastric Cancer" ~ 0,
    primary_disease == "Gallbladder Cancer" ~ 0,
    primary_disease == "Esophageal Cancer" ~ 0,
    primary_disease == "Colon/Colorectal Cancer" ~ 0,
    primary_disease == "Liver Cancer"~ 0,
    
    primary_disease == "Ovarian Cancer" ~ 1,
    primary_disease == "Cervical Cancer" ~ 1,
    primary_disease == "Endometrial/Uterine Cancer" ~ 1,
    primary_disease == "Teratoma" ~ 1,
    
    primary_disease == "Bone Cancer" ~ 2,
    primary_disease == "Skin Cancer" ~ 2,
    primary_disease == "Fibroblast" ~ 2,
    primary_disease == "Sarcoma" ~ 2,
    primary_disease == "Liposarcoma" ~ 2,
    
    primary_disease == "Brain Cancer" ~ 3,
    primary_disease == "Neuroblastoma" ~ 3,
    primary_disease == "Rhabdoid" ~ 3,
    
    primary_disease == "Breast Cancer" ~ 4,
    
    primary_disease == "Head And Neck Cancer" ~ 5,
    primary_disease == "Thyroid Cancer" ~ 5,
    
    primary_disease == "Leukemia"~ 6,
    primary_disease == "Lymphoma"~ 6,
    primary_disease == "Myeloma"~ 6,
    
    primary_disease == "Bladder Cancer" ~ 7,
    primary_disease == "Kidney Cancer" ~ 7,
    primary_disease == "Prostate Cancer" ~ 7,
    
    primary_disease == "Lung Cancer" ~ 8,
    
    TRUE ~ 9 
  ))


label = trunc_df2_bis$CancerType
df1_bis = df1 %>%
  filter(DepMap_ID %in% trunc_df2_bis$DepMap_ID)
multiclass_dataset = cbind(df1_bis, label)

# training and test sets
set.seed(8675309)
n.train = floor(.80*dim(df1_bis)[1])
training = sample(1:dim(df1_bis)[1], size = n.train, replace = FALSE)

multiclass.training = multiclass_dataset[training, ]  
multiclass.test = multiclass_dataset[-training, ]


# write.csv(multiclass.training, "dataset/multiclass_training_no_weird_obs.csv", rownames = F)
# write.csv(multiclass.test, "dataset/multiclass_test_no_weird_obs.csv", rownames = F)


## multiclass
trunc_df2_bis = trunc_df2 %>%
  filter(primary_disease != "Non-Cancerous") %>%
  mutate(CancerType = case_when (
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
    
    primary_disease == "Eye Cancer" ~ 0,
    (primary_disease == "Engineered") & (sample_collection_site == "Eye") ~ 0,
    
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

# training and test sets
# write.csv(multiclass.training, "dataset/multiclass_training.csv", rownames = F)
# write.csv(multiclass.test, "dataset/multiclass_test.csv", rownames = F)
# write.csv(multiclass_dataset,"dataset/multiclass_dataset.csv", rownames = F)


# plot of the 10 classes
ggplot(trunc_df2_bis, aes(x = as.factor(CancerType))) +
  geom_bar(aes(fill = primary_disease), width = 0.55, 
           position = position_dodge(width = 0.75)) +
  ggtitle("Cancers by Body Location/System") +
  labs(x = "Class", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Eye", "Gastrointestinal", "Gynecologic",
                              "Musculoskeletal", "Neurologic", "Breast",
                              "Head-Neck", "Hematologic", "Genitourinary", "Lung")) +
  scale_fill_discrete(name = "Primary disease")


p = ggplot(trunc_df2_bis, aes(x = as.factor(CancerType))) +
  geom_bar(aes(fill = as.factor(CancerType)), width = 0.55, show.legend = FALSE) +
  ggtitle("Cancers by Body Location/System") +
  labs(x = "Class", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Eye", "Gastrointestinal", "Gynecologic",
                              "Musculoskeletal", "Neurologic", "Breast",
                              "Head-Neck", "Blood", "Genitourinary", "Lung")) +
  scale_fill_discrete(name = "Primary disease") +
  theme(panel.border = element_blank(), #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))
p

# save_plot("report/plot1.png", p)  

library(stringi)
library(ggplot2)
library(dplyr)
library(knitr)
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

ggplot(nas, aes(x = row)) +
  geom_bar(stat = "count") +
  ggtitle("How many NA per row")

trunc_df2 %>%
  filter(DepMap_ID %in% row_na) %>%
  select(DepMap_ID, primary_disease, sample_collection_site, lineage) %>%
  kable()
# Nothing particular, obs seem unrelated 

NaCount = nas %>%
  group_by(row) %>%
  summarise(count = n())
NaCount = cbind(df1$DepMap_ID[unique(nas$row)], NaCount)
colnames(NaCount) = c("DepMap_ID", "Row_index", "Count")
kable(NaCount)

NaCount[NaCount$Count > 1000, ]$Row_index # to be deleted

new_df1 = df1 %>%
  filter(DepMap_ID %notin% NaCount[NaCount$Count > 1000, ]$DepMap_ID)

# ggplot(data.frame(which(is.na(new_df1), arr.ind = TRUE)), 
#        aes(x = row)) +
#   geom_bar(stat = "count") +
#   ggtitle("How many NA per row")

# And let us replace NA with 0
new_df1[is.na(new_df1)] = 0


## check of row sum
sum(new_df1[1, -1])
sum(new_df1[1000, -1])
sum(new_df1[504, -1])
sum(new_df1[750, -1])
sum(new_df1[9, -1]) #
sum(new_df1[379, -1]) #
sum(new_df1[417, -1]) #
sum(new_df1[548, -1]) #
sum(new_df1[582, -1]) #


# check of col sum
sum(new_df1[, 2])
sum(new_df1[, 121])
sum(new_df1[, 16897])
sum(new_df1[, 5598])
sum(new_df1[, 489])
sum(new_df1[, 4970])
sum(new_df1[, 24])
sum(new_df1[, 697])


# more_than_one = select_if(new_df1, ~any(. > 1.01))

x = seq(0, 1027)
for (i in 1:1027) {
  x[i] = length(which(new_df1[i, -1] > 1.0))
}
sum(x) # 1027
x[1028] # = 1027

############################################################
# ♫ # Comparison btw primary_disease and sample_collection_site
############################################################

cat("n. diseases =", length(unique(trunc_df2$primary_disease)), 
    "VS. n. sites =", length(unique(trunc_df2$sample_collection_site)))
# more sites then disease: let us investigate more and see if there is 
# an evident substructure

CancerSiteCount = trunc_df2 %>%
  select(DepMap_ID, primary_disease, sample_collection_site) %>%
  rename(disease = primary_disease,
         site = sample_collection_site)

CancerSiteCount %>%
  group_by(disease, site) %>% 
  summarise(count = n()) %>%
  kable()
# Stranger things: probably it is not a good idea looking at the sample site for the clusters

ggplot(CancerSiteCount, aes(x = disease)) +
  geom_histogram(aes(fill = site), stat = "count") +
  ggtitle("Number of cancer per site") +
  labs(x = "Type of cancer", y = "How many obs") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        plot.title = element_text(hjust = 0.5))
# CONCLUSION: just stick to the primary_disease for clusters 


############################################################
# ♫ # Count obs wrt primary disease
############################################################
trunc_df2 = df2 %>%
  filter(DepMap_ID %in% df1$DepMap_ID)

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

# ♫ # Gastroinstestinal cancers VS. All the others
trunc_df2 = df2 %>%
  filter(DepMap_ID %in% df1$DepMap_ID)

trunc_df2 = trunc_df2 %>%
  mutate(isGastro = case_when (
    primary_disease ==  "Bile Duct Cancer" ~ 1,
    primary_disease == "Kidney Cancer"~ 1,
    primary_disease == "Gastric Cancer"~ 1,
    primary_disease == "Gallbladder Cancer"~ 1,
    primary_disease == "Esophageal Cancer"~ 1,
    primary_disease == "Colon/Colorectal Cancer"~ 1,
    primary_disease == "Liver Cancer"~ 1,
    TRUE~ 0 
  ))

# appending this column to the CRISP dataset to have labeled data
label = trunc_df2$isGastro
gastro_dataset = cbind(df1, label)

# and split the data into training and test 
set.seed(8675309)
n.train = floor(.80*1032)
training = sample(1:1032, size = n.train, replace = FALSE)

gastro.training = gastro_dataset[training, ]  
gastro.test = gastro_dataset[-training, ]

set.seed(8675309)


# double check:
#df2 %>%
#  filter(DepMap_ID %in% df1$DepMap_ID, primary_disease == "Lung Cancer") %>%
#  dim()
# same


############################################################
# ♫ # Lung cancer VS. All
############################################################

# Adding a column containing whether the obs is a Lung cancer or not
#df2 = df2 %>%
# mutate(isLung = as.numeric(primary_disease == "Lung Cancer"))

# appending this column to the CRISP dataset to have labeled data
#label = df2[rownames(df1),]$isLung
#lung_dataset = cbind(df1, label)

# and split the data into training and test 
#set.seed(8675309)
#n.train = floor(.80*1032)
#training = sample(1:1032, size = n.train, replace = FALSE)

#lung.training = lung_dataset[training, ]  
#lung.test = lung_dataset[-training, ]

# RMK: Lung cancer is the 12.3% of the total obs --> ? minority class ?

# write.csv(lung.training,"~/Documents/sds/sml/SML-project/dataset/lung_training.csv", 
#          row.names = TRUE)
# write.csv(lung.test,"~/Documents/sds/sml/SML-project/dataset/lung_test.csv", 
#          row.names = TRUE)

#####
# ♫ # Without weird obs:
#####

weird_obs = c(which(trunc_df2$primary_disease == "Non-Cancerous"), 
              which(trunc_df2$primary_disease == "Engineered"))
# trunc_df2_bis = trunc_df2[-weird_obs, ]
df1_bis = df1[-weird_obs, ]
rownames(df1_bis) = df1_bis$DepMap_ID

# Lung cancer VS. Al
#label_bis = df2[rownames(df1_bis),]$isLung
#lung_dataset_bis = cbind(df1_bis, label_bis)

#set.seed(8675309)
#training_bis = sample(1:1024, size = floor(.80*1024), replace = FALSE)
#lung_bis.training = lung_dataset_bis[training_bis, ]  
#lung_bis.test = lung_dataset_bis[-training_bis, ]


############################################################
# ♫ # 3/4 classes processing
############################################################

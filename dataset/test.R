library(stringi)
library(ggplot2)
library(dplyr)
library(knitr)
#####################################################

df1 = read.csv("CRISPR_gene_dependency.csv")
df2 = read.csv("sample_info.csv")

df1$DepMap_ID %in% df2$DepMap_ID
length(which(df1$DepMap_ID %in% df2$DepMap_ID)) # as the n of df1

# ♠ RMK: lines from 17 up to 64 deal with the whole df2, so that we have a
# ♠      general overview


# ♫ # (de-)Capitalize columns
unique(df2$primary_disease)
df2$primary_disease = stri_trans_totitle(df2$primary_disease)
unique(df2$primary_disease)

unique(df2$sample_collection_site)
df2$sample_collection_site = chartr(old = "_", new = " ", df2$sample_collection_site)
df2$sample_collection_site = stri_trans_totitle(df2$sample_collection_site)
unique(df2$sample_collection_site)


# ♫ # Looking for weird obs labels
length(which(df2$primary_disease == "Unknown"))
which(df2$primary_disease == "Unknown") # obs >1032 --> not in df1

length(which(df2$primary_disease == ""))
which(df2$primary_disease == "") # obs >1032 --> not in df1

length(which(df2$primary_disease == "Non-Cancerous"))
which(df2$primary_disease == "Non-Cancerous") # obs >1032 --> not in df1

length(which(df2$primary_disease == "Immortalized"))
which(df2$primary_disease == "Immortalized") # obs >1032 --> not in df1

length(which(df2$primary_disease == "Engineered"))
which(df2$primary_disease == "Engineered") # 5 obs are in df1
df2$sample_collection_site[c(49, 64, 170, 494, 641)]
# just 5 and each from different site --> ok to delete


# RMK: also lineage can be used as Cancer type classifications, but seems 
#      more complicated (since more details)
# --> idea: use it as double check for engineered:
df2$lineage = chartr(old = "_", new = " ", df2$lineage)
df2$lineage = stri_trans_totitle(df2$lineage)
unique(df2$lineage)
# More obs than primary_disease: just with a quick look, we see engineered
# is divided in about 10 other categories 
which(startsWith(df2$lineage, "Engineered"))
# Same 5 initial obs as before --> ok!
which(startsWith(df2$lineage, "Engineered")) %in% which(df2$primary_disease == "Engineered")
# only last obs do not corresponds to primary disease
# CONLUSION: it is ok to delete those engineered obs (if it is the case)


# ♫ # Comparison btw primary_disease and sample_collection_site
cat("n. diseases =", length(unique(df2$primary_disease)), 
    "VS. n. sites =", length(unique(df2$sample_collection_site)))

CancerSiteCount = df2[1:1032, ] %>%
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


# ♫ # Count obs wrt primary disease
trunc_df2 = df2[1:1032,]
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


###### Update 08/11/2021
# ♫ # Lung cancer VS. All

# Add a column containing whether the obs is a Lung cancer or not
df2$primary_disease = stri_trans_totitle(df2$primary_disease)
lung_df2 = df2[1:1032, ] %>%
  mutate(isLung = (primary_disease == "Lung Cancer"))

# and split the data into training and test 
set.seed(8675309)
n.train = floor(.80*1032)
training = sample(1:1032, size = n.train, replace = FALSE)

lung.training = lung_df2[training, ]  
lung.test = lung_df2[-training, ]

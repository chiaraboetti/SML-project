df1 = read.csv("CRISPR_gene_dependency.csv")
df2 = read.csv("sample_info.csv")


df1$DepMap_ID %in% df2$DepMap_ID
length(which(df1$DepMap_ID %in% df2$DepMap_ID)) # as the n of df1

# (de-)Capitalize columns
unique(df2$primary_disease)
unique(df2$sample_collection_site)

library(stringi)
df2$primary_disease = stri_trans_totitle(df2$primary_disease)
unique(df2$primary_disease)
df2$sample_collection_site = stri_trans_tolower(df2$sample_collection_site)
unique(df2$sample_collection_site)


# Looking for weird obs labels
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


# RMK: also lineage can be used as Cancer type classifications, but seems more complicated
# --> idea: use it as double check (to do)

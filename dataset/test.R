df2 = read.csv("CRISPR_gene_effect.csv")
df3 = read.csv("sample_info.csv")


df2$DepMap_ID %in% df3$DepMap_ID
unique(df3$primary_disease)


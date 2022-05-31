require(DESeq2)
require(tidyverse)
library(ggplot2)
require(httr)
require(jsonlite)
library(plyr)
df <- read.csv("../data/data_raw_counts.csv",header = TRUE, sep = ",",row.names=NULL)

###
# Multiple IDs to convert - use a POST request
###
url = "https://biotools.fr/human/ensembl_symbol_converter/"
ids = c(df$ensgene)
ids_json <- toJSON(ids)

body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)

output = fromJSON( content(r, "text"), flatten=TRUE)
gene_id <- ldply(output)
gene_id
colnames(gene_id) <- c("row","gene")
metadata <- read.csv("../data/covariates.csv")
case_control<-metadata[,c("Participant_ID","Case_Control")]
dds <- DESeqDataSetFromMatrix(countData=df, 
                              colData=case_control, 
                              design=~Case_Control, tidy = TRUE)
dds
dds <- DESeq(dds)
res <- results(dds, tidy=TRUE)
head(res)
summary(res)

## Replace ENS with Gene ID 
res<-merge(res,gene_id,by="row")
write.csv(res,"../data/data_de.csv", row.names = FALSE)

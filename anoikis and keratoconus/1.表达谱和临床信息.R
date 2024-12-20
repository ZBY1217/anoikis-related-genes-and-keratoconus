# -------------------导入函数-------------------
# source("src/0.config.R")
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources, source))
library(tidyverse)
library(data.table)
# -------------------设置结果路径-------------------
od <- "data"
# -------------------------------获取表达信息和临床信息-------------------------------
load("/Pub/Users/cuiye/database/gene_map.RData")
save(gene_map, file = str_glue("{od}/gene_map.RData"))



#gse204839
expression = fread(file = "data/GSE204839GSE204839_GPL21185/raw_expression.tsv")

expression = expression %>% column_to_rownames("Gene_Symbol")
save(expression, file = "data/expression.RData")

clinical = fread(file = "data/GSE204839GSE204839_GPL21185/raw_clinical.tsv")

clinical = clinical %>% select(c("geo_accession", "age:ch1", "disease:ch1", "gender:ch1"))

save(clinical, file = "data/clinical.RData")


#GSE77938
expression_1 = fread(file = "data/GSE77938/GSE77938_discovery_gene_tpm.txt")
expression_2 = fread(file = "data/GSE77938/GSE77938_replication_gene_tpm.txt")
expression_all = inner_join(expression_1, expression_2)
expression_all = expression_all %>% select(!c(1,3)) %>% group_by(`Gene symbol`) %>% 
             summarise_all(median, na.rm=TRUE) %>% column_to_rownames("Gene symbol")
save(expression_all, file = "data/expression_GSE77938_all.RData")
expression = expression %>% filter(`Gene biotype` == "protein_coding") %>% 
             select(!c(1,3)) %>% group_by(`Gene symbol`) %>% 
             summarise_all(median, na.rm=TRUE) %>% column_to_rownames("Gene symbol")
save(expression, file = "data/expression_GSE77938.RData")

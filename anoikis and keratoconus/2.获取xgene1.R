# -------------------导入函数-------------------
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources, source))
library(tidyverse)
library(data.table)
# -------------------设置结果路径-------------------
od <- "data/"

# -------------------------------获取xgene-------------------------------
# *******目的：从下载的基因集中，得到基因的symbol向量（存为xgene）*******
source("/Pub/Users/cuiye/RCodes/UserCode/newlover/key_match_path_ingmt.R")
source("/Pub/Users/cuiye/RCodes/UserCode/newlover/path_match_gene_ingmt.R")
path <- key_match_path_ingmt(key = "ANOIKIS", database = "GOBP")
xgene <- path_match_gene_ingmt(path = path[1])

# ---------------------------xgene与表达谱取交集------------------------------
load("/Pub/Users/renxz/project/F221009004/data/expression_GSE77938.RData")
length(xgene)
xgene <- intersect(rownames(expression), xgene)
length(xgene)
save(xgene, file = str_glue("{od}xgene1.RData"))
load("data/xgene.RData")


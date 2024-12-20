library(data.table)
library(tidyverse)
# 设置结果路径
od = "data"

# -------------------------------获取xgene-------------------------------
# *******目的：从下载的基因集中，得到基因的symbol向量（存为xgene）*******


# 以下是参考代码 分别是https://www.gsea-msigdb.org/gsea/index.jsp
# 和http://amigo.geneontology.org/amigo 数据库下载的数据处理
# 1.msigdb
p_load(rjson)
json = fromJSON(file = "data/REACTOME_GLUCOSE_METABOLISM.v7.5.1.json")
xgene = json$REACTOME_GLUCOSE_METABOLISM$geneSymbols
save(xgene,file="data/xgene.RData")
#2.amigo 
xgene = fread("data/xgene.csv",sep=':',header = F) %>% pull(V2)
suppressPackageStartupMessages(library(clusterProfiler))
conflict_prefer("first", "dplyr")
xgene = bitr(geneID=xgene,fromType="UNIPROT",toType="SYMBOL",OrgDb="org.Hs.eg.db")
xgene %<>% pull(2)
length(xgene)
#tab sep file
xgene <- fread("data/osr_77.txt", skip = 1) %>% pull(1)

# ---------------------------xgene与表达谱取交集------------------------------
load("data/train_data.RData")
length(xgene)

xgene <- intersect(rownames(train_data$data_exprs),xgene)
length(xgene)
save(xgene,file = str_glue("{od}/xgene.RData"))



load("data/xgene.RData")
xgene

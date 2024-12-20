rm(list = ls())
library(pacman)
od = "data/GSE204839"
dir.create("data/GSE204839", recursive = T)
GEO <- "GSE204839"
p_load(data.table, GEOquery, tidyverse, limma, AnnoProbe, ggsignif, statmod, parallel,magrittr)
#--------------------------------1.下载-----------------------------------------
gset <- getGEO(GEO, destdir = "data/GSE204839", getGPL = T)
#--------------------------------2.整理-----------------------------------------
# extract the expression matrix and phenotype data
eSet <- gset[[2]] # GPL570
plt <- eSet$platform_id[1]
probes_expr <- exprs(eSet)
dim(probes_expr)
head(probes_expr[, 1:4])
range(probes_expr, na.rm = TRUE)
probes_exprs <- probes_expr %>% as.data.table(keep.rownames = "probe_id")
probes_exprs[1:3, 1:3]
#----------------------------3.获取探针注释文件---------------------------------
idmaps <- function(ann_file, ID = "ID", Symbol = "Gene symbol", skip = 27, sep = "\t") {
    temp <- fread(ann_file, sep = sep, skip = skip)
    vars <- c(ID, Symbol)
    temp <- temp[, ..vars]
    setnames(temp, c("probe_id", "symbol"))
    temp %<>% separate_rows(symbol, sep = "///")
    temp <- temp[!is.null(temp$symbol), ]
    temp <- temp[!is.na(temp$symbol), ]
    temp <- temp[temp$symbol != "", ]
    temp$symbol %<>% str_remove_all(" ") %>% str_remove_all("\"")
    return(temp)
}
probe2symbol <- idmaps("data/GSE204839/GPL21185-21174.txt", Symbol = "GENE_SYMBOL", skip = 16)
probe2symbol[1:10, 1:2]
#----------------------------4.ID转换-------------------------------------------
transid <- function(probe2symbol, exprSet, method = "median") {
    probe2symbol$probe_id %<>% as.character()
    exprSet$probe_id %<>% as.character()
    ex <- probes_exprs %>%
        dplyr::inner_join(probe2symbol, by = "probe_id") %>% # 合并探针的信息
        dplyr::select(-probe_id) %>% # 去掉多余信息
        dplyr::select(symbol, everything()) %>% # 重新排列，
        dplyr::mutate(ref = apply(across(where(is.numeric)), 1, method)) %>%
        dplyr::arrange(desc(ref)) %>% # 把表达量的平均值按从大到小排序
        dplyr::select(-ref) %>% # 反向选择去除rowMean这一列
        dplyr::distinct(symbol, .keep_all = T)
    return(ex)
}

genes_expr <- transid(probe2symbol, probes_exprs, "mean") %>%
    tidyr::drop_na() %>%
    column_to_rownames(var = "symbol")

genes_expr[1:10, 1:3]

dir.create(str_glue("{od}{GEO}_{plt}/"), recursive = T)
fwrite(as.data.table(genes_expr,keep.rownames = "Gene_Symbol"),file=str_glue("{od}{GEO}_{plt}/raw_expression.tsv"),sep='\t')

# clinical
fwrite(pData(eSet),file = str_glue("{od}{GEO}_{plt}/raw_clinical.tsv"),sep='\t')

load("/Pub/Users/wangyk/project/Poroject/F221009004/data/xgene.RData")

xgene %>% paste0(collapse = '、') %>% cat('\n')


df <- read.delim("/Pub/Users/wangyk/project/Poroject/F221009004/results/2.WGCNA/res/merged_infor.txt")

dim(df)
head(df)

rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))
library(tidyverse)
library(magrittr)
conflict_prefer("between", "dplyr",quiet = TRUE)

# 设置结果路径
od <- "results/1.cluster/"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/expression_GSE77938.RData")
load("results/1.cluster/cluster.RData")
expression = log(expression + 1)
dat = expression[,rownames(cluster_hc_maximum_2)]
cluster_hc_maximum_2$cluster = str_c("cluster", cluster_hc_maximum_2$cluster)
# 免疫细胞浸润

immunescore <- immune_score(exp = dat, arrays = FALSE, perm = 200, group_list = NULL,
                         tumor = TRUE, scale_mrna = TRUE, od = od, ssgsea_geneSets = NULL,
                         method = c("cibersort", "xcell", "estimate", "ssgsea"))
save(immunescore, file = str_glue("{od}/immunescore.RData"))


# Adaptive
type28 <- openxlsx::read.xlsx("/Pub/Users/cuiye/database/geneSet/28type_from_dio10.1016_j.celrep.2016.12.019_for_ssGSEA.xlsx", sheet = 1, colNames = T)
data <- type28 %>% 
    filter(Immunity == "Adaptive") %>%
    pull(Cell.type) %>%
    unique()
dir.create(str_glue("{od}immune_infiltration"))
characteristics_plot_by_group(
    characteristics_score = immunescore$ssgsea %>% select(sample, data),
    Group = cluster_hc_maximum_2,
    od = str_glue("{od}immune_infiltration/ssgsea_adaptive"),
    type = "ssgsea",
    heatplot_by_scale = TRUE,
    cluster_rows = TRUE
)

# Innate
data <- type28 %>% 
    filter(Immunity == "Innate") %>%
    pull(Cell.type) %>%
    unique()
characteristics_plot_by_group(
    characteristics_score = immunescore$ssgsea %>% select(sample, data),
    Group = cluster_hc_maximum_2,
    od = str_glue("{od}immune_infiltration/ssgsea_innate/"),
    type = "ssgsea",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

# other
lapply(c("cibersort", "xcell", "estimate"), function(x) {
    characteristics_plot_by_group(
        characteristics_score = immunescore[[x]],
        Group = cluster_hc_maximum_2,
        od = str_glue("{od}immune_infiltration/other/"), type = x,
        heatplot_by_scale = TRUE,
        cluster_rows = TRUE
    )
})

#HLA差异
HLA_gene = rownames(expression)[rownames(expression) %>% str_detect("^HLA")]

source("src/functions/v_characteristics_plot_by_group.R")
v_characteristics_plot_by_group(
    characteristics_score = dat[HLA_gene, ] %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = cluster_hc_maximum_2 %>% rename(Group = cluster),
    od = str_glue("{od}immune_infiltration/HLA/"),
    type = "ssgsea",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

#趋化因子
chemokine = read.table("data/趋化因子.txt") %>% pull(1)
chemokine = intersect(chemokine, rownames(expression))

v_characteristics_plot_by_group(
    characteristics_score = dat[chemokine, ] %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = cluster_hc_maximum_2 %>% rename(Group = cluster),
    od = str_glue("{od}immune_infiltration/chemokine/"),
    type = "ssgsea",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

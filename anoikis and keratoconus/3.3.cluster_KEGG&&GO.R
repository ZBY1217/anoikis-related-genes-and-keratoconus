rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))
library(clusterProfiler)
library(GSVA)
library(data.table)
#结果路径
od = "results/1.cluster/kegg_GO"
dir.create(od, recursive = T)

#数据处理
load("data/xgene1.RData")
load("results/1.cluster/cluster.RData")
load("data/expression_GSE77938.RData")

##表达数据
expression = log(expression + 1)
dat = expression[,rownames(cluster_hc_maximum_2)]

##通路基因
geneSets <- GSEABase::getGmt("/Pub/Users/wangcy/database/MsigDB/msigdb.v2022.1.Hs.symbols.gmt")
aa = map_dfr(seq_len(length(geneSets)), function(x){
    path_name = geneSets[[x]]@setName
    path_name = data.frame(num = x, pathname = path_name)
})
aa = aa %>% filter(str_detect(aa$pathname, "(KEGG|GOBP)"))
geneSets = geneSets[c(aa$num)]

#计算ssgsea
gsva_matrix<- gsva(as.matrix(dat), geneSets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

#绘图
source("src/functions/Heatmap_manul.R",local = TRUE)
group_infor = cluster_hc_maximum_2 %>% rename("Cluster" = "cluster")
set.seed(123456)


###kegg
kegg = as.data.frame(gsva_matrix)[str_detect(rownames(gsva_matrix), "KEGG"), ]

Heatmap_manul(
    data_input = kegg,
    color_used = color_fun1,
    group_infor = group_infor,
    Colored_OtherInfor = color_fun1,
    show_rownames = T,
    saveplot = TRUE, 
    output_dir = od,
    var_name = "KEGG", 
    width_used = 12, 
    height_used = 10
)

result_kegg = characteristics_plot_by_group(
    characteristics_score = kegg %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = group_infor,
    od = od, type = "ssgsea",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
source("src/functions/v_characteristics_plot_by_group.R")
path = result_kegg %>% filter(p.signif !=" ") %>% rownames()
result_kegg1 = v_characteristics_plot_by_group(
    characteristics_score = kegg %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = group_infor,
    od = od, type = "kegg",
    heatplot_by_scale = TRUE, cluster_rows = TRUE,
    feature2show = path
)

###GO
go = (gsva_matrix %>% as.data.frame())[str_detect(rownames(gsva_matrix), "GO"), ]
result_go = characteristics_plot_by_group(
    characteristics_score = go %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = group_infor,
    od = od, type = "go",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
# source("src/functions/v_characteristics_plot_by_group.R")
path = fread("results/1.cluster/kegg_GO/SupplementaryTable_go_sign_stat.txt")
path = (path %>% filter(p.signif !=" ") %>% arrange(pvalue))[1:20, ]
result_go1 = v_characteristics_plot_by_group(
    characteristics_score = go[path$V1, ] %>% t() %>% data.frame() %>% rownames_to_column("sample"),
    Group = group_infor,
    od = od, type = "go",
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

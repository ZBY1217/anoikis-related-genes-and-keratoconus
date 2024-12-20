rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))

#结果路径
od = "results/1.cluster/pca"
dir.create(od, recursive = T)

#数据处理
load("data/xgene1.RData")
load("results/1.cluster/cluster.RData",verbose = T)

# load("results/1.cluster/cluster_hc_spearman_3.RData")
load("data/expression_GSE77938.RData")
expression = log(expression + 1)

#pca
dat = expression[xgene,rownames(cluster_hc_maximum_2)] %>% na.omit() %>% t()
p_load(FactoMineR,factoextra)
pca = PCA(dat, scale.unit = TRUE, ncp = 2, graph = F)
p = fviz_pca_ind(X = pca, ## pca对象
  axes = 1:2, ## 展示的两个主成分
  geom = 'point', ## 展示individual的形式
  habillage = factor(cluster_hc_maximum_2$cluster), ## individual用来分组的变量
  legend.title = 'Groups', ## 分组变量的title
  palette = 'lancet', ## 颜色面板
  addEllipses = T,  ## 是否绘制椭圆
  ellipse.level = 0.7, ## 椭圆的大小
  title = 'PCA Plot', ## 标题
  mean.point = F ## 不删除每个组的重心
  )+
  theme(plot.title = element_text(hjust = 0.5,size=20))
ggsave(plot = p,str_c(od,"xgene_PCA2d_cluster.pdf"),width=4.7,height=4.7)

gr <- factor(xgene_cluster$Cluster) 
save(pca,gr,file = str_glue("{od}cluster_pca.RData"))


#xgene在不同分组中表达情况
#结果路径
od = "results/1.cluster/xgene"
dir.create(od, recursive = T)
#分组信息
#热图
# group_infor <- merge(train_data$data_clinical, xgene_cluster) %>%
#     column_to_rownames("sample") %>%
#     dplyr::select(Cluster, any_of(parameter_list$clinical_feature_list))
group_infor = cluster_hc_maximum_2 %>% rename("Cluster" = "cluster")

source("src/functions/Heatmap_manul.R",local = TRUE)

set.seed(123456)
Heatmap_manul(
    data_input = expression[xgene, ] %>% na.omit(),
    color_used = color_fun1,
    group_infor = group_infor,
    Colored_OtherInfor = color_fun1,
    show_rownames = ifelse(length(xgene)>50,FALSE,TRUE),
    saveplot = TRUE, 
    output_dir = od,
    var_name = "xgene", 
    width_used = 12, 
    height_used = 8
)

#箱线图
characteristics_score <- expression[xgene, ] %>%
    na.omit() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")

# 不同临床特征分组

gr <- group_infor %>%
    dplyr::rename(Group = Cluster) %>%
    mutate(Group = factor(Group))
source("src/functions/v_characteristics_plot_by_group.R",local = TRUE)
v_characteristics_plot_by_group(
    characteristics_score = characteristics_score,
    feature2show=33,
    Group = gr,
    od = od # str_glue("{od}/{x}/"),
    , type = "Cluster"
)

#deg
od = "results/1.cluster/deg"

#数据处理
expression = log(expression+1)
dat = expression[,rownames(cluster_hc_maximum_2)]
sample_pdata = data.frame(Group = str_c("Cluster", cluster_hc_maximum_2$cluster)) 
rownames(sample_pdata) = rownames(cluster_hc_maximum_2)
sample_pdata = sample_pdata %>% rownames_to_column("sample")

#差异分析
deg_res <- limma_deg(
    od = od, DEG_exp = dat, DEG_pdata = sample_pdata,
    controlLabel = "Cluster1", caseLabel = "Cluster2",
    DEG_FC = 1, DEG_P = 0.05, pvalue = NULL, saveplot = T, color_fun = c("#E41A1C","#377EB8" )
)

save(deg_res, file = str_glue("{od}/deg.RData"))


###多分型
# unique_cluster <- unique(sample_pdata$Group) %>% sort()
# combine_deg_res <- map_dfr(unique_cluster, function(cluster) {
#     sample_pdata <- sample_pdata %>% dplyr::mutate(Group = ifelse(Group == cluster, cluster, "Rest"))
#     deg_res <- limma_deg(
#         od = od, DEG_exp = dat, DEG_pdata = sample_pdata,
#         controlLabel = "Rest", caseLabel = cluster,
#         DEG_FC = 0.585, DEG_P = 0.05, pvalue = NULL, saveplot = T, color_fun = color_fun1
#     )
#     return(deg_res$nrDEG %>% rownames_to_column("gene") %>% dplyr::mutate(Cluster = cluster))
# })


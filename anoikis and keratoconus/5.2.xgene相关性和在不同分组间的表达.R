rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))
# source("src/0.config.R")
# conflict_prefer("between", "dplyr",quiet = TRUE)
# 导出路径
od <- "results/3.xgene_landscape/exprs/"
suppressWarnings(dir.create(od,recursive=TRUE))
library(pacman)
p_load("ggplotify", "tidyverse", "magrittr")
# 导入数据
load("data/expression_GSE77938.RData")
load("data/xgene1.RData")
expression = log2(expression+1)
characteristics_score <- expression[xgene, ] %>%
    na.omit() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")

pdata = data.frame(sample = characteristics_score$sample,
                   Group = ifelse(str_detect(characteristics_score$sample, "KC"), "disease", "healthy"))


# 正常vs对照
source("src/functions/v_characteristics_plot_by_group.R", local = TRUE)
if (all(table(pdata$Group) > 5) && unique(pdata$Group)>1) {
    v_characteristics_plot_by_group(
        characteristics_score = characteristics_score, feature2show = 40,
        Group = pdata %>% column_to_rownames("sample"), od = str_glue("{od}/Tumor_vs_Normal/"), type = "xgene"
    )
}

#差异分析
od = str_glue("{od}/Tumor_vs_Normal/")
dat = expression[xgene, ]
deg_res <- limma_deg(
    od = od, DEG_exp = dat, DEG_pdata = pdata,
    controlLabel = "healthy", caseLabel = "disease",
    DEG_FC = 0.585, DEG_P = 0.05, pvalue = NULL, saveplot = T, color_fun = c("#E41A1C","#377EB8" )
)





#相关性
##加载包
# devtools::install_github("Hy4m/linkET", force = TRUE)
# library(linkET)
library(corrplot)
##加载数据
load("data/expression_GSE77938.RData")
load("data/xgene1.RData")
expression = log2(expression+1)
expression_xgene = expression[xgene, ]
expression_disease_xgene = expression_disease[xgene, ]
expression_disease = expression %>% select(starts_with("KC"))

##生成颜色

col3 = colorRampPalette(c("blue", "white", "red"))
col4 = colorRampPalette(c("green", "white", "purple"))

##显著性检验
all_p = cor.mtest(t(expression_xgene), conf.level = .95)
disease_p = cor.mtest(t(expression_disease_xgene), conf.level = .95)

##绘图
pdf(file = "results/1.xgene_landscape/cor.pdf")
corrplot(cor(t(expression_xgene)), type='lower', method='circle', tl.pos='lt', cl.cex=1,tl.cex = 1,
    tl.col = "black",col = col3(100), diag= T, mar = c(2, 2, 2, 4), p.mat = all_p$p, sig.level = 0.05,
    pch.cex = 0.5)
corrplot(cor(t(expression_disease_xgene)),add=TRUE, type='upper', method='square', tl.pos='n',
    diag=F, cl.cex=1, col = col4(100), p.mat = disease_p$p, sig.level = 0.05, pch.cex = 0.5)
dev.off()



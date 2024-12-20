rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))
library(tidyverse)
#结果路径
od = "results/2.WGCNA/"
dir.create(od, recursive = T)

#数据处理
load("data/xgene1.RData")
load("results/1.cluster/cluster.RData")
load("data/expression_GSE77938.RData")
expression = log2(expression + 1)
dat = expression[,rownames(cluster_hc_maximum_2)]

sample_data = data.frame(sample = rownames(cluster_hc_maximum_2), 
                        cluster1 = ifelse(cluster_hc_maximum_2$cluster == 1, 1, 0),
                        cluster2 = ifelse(cluster_hc_maximum_2$cluster == 2, 1, 0))

powers = c(c(1:10), seq(from = 12, to=20, by=2))
datExpr = t(dat)


library(WGCNA)

gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers1 <- c(seq(1, 10, by=1), seq(12, 20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers1)
RpowerTable <- pickSoftThreshold(datExpr, powerVector = powers1)[[2]]

cex1 = 0.7
par(mfrow = c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], xlab = "soft threshold (power)", ylab = "scale free topology model fit, signes R^2", type = "n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels = powers1, cex = cex1, col = "red")
abline(h = 0.95, col = "red")
plot(RpowerTable[,1], RpowerTable[,5], xlab = "soft threshold (power)", ylab = "mean connectivity", type = "n")
text(RpowerTable[,1], RpowerTable[,5], labels = powers1, cex = cex1, col = "red")

sft$powerEstimate






#
qc_res <- wgcna_qc(
    exp = dat,
    pheno = sample_data,
    method = "average", cutHeight = 150, cutMad = 0, od = paste0(od, "QC"), width = 9, height = 6
)
allowWGCNAThreads(nThreads = 20)
enableWGCNAThreads(nThreads = 20)
WGCNAnThreads()

softpower_res <- wgcna_picksoftpower(exp = qc_res$use_exp, od = paste0(od, "res"), net_type = "unsigned", Rsquared_cut = 0.85, width = 8, height = 6)


wgcna_net <- wgcna_built_net(exp = qc_res$use_exp, pheno = qc_res$use_pheno, power = 6, od = paste0(od, "res"), height = 15, width = 6)
write.table(wgcna_net$merged_infor, file = str_glue("{od}/res/merged_infor.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)


# 根据选定参数筛选WGCNA网络模块基因（核心基因），生成Figure_{select_pheno}_{select_module}_cor.pdf
# 读取模块数
colorcount <- read.delim(paste0(od, "res/SupplementaryTable_module_gene_count.txt")) 

# 读取相关性值，并加入其模块数
cor <- read.delim(paste0(od, "res/SupplementaryTable_modTraitCor.txt")) %>%
    rename("corvalue" = "cluster2") %>%
    rownames_to_column("module") %>%
    mutate(module = str_replace(module, "ME", "")) %>%
    inner_join(colorcount, by = "module")

# 读取p值，并加入其模块数
cor_p_count <- read.delim(paste0(od, "res/SupplementaryTable_modTraitP.txt")) %>%
    rename("pvalue" = "cluster2") %>%
    rownames_to_column("module") %>%
    select(c("module", "pvalue")) %>%
    mutate(module = str_replace(module, "ME", "")) %>%
    inner_join(cor, by = "module") %>%
    rename("count" = "number") %>%
    select(c("module", "corvalue", "pvalue", "count", everything()))

# view(cormax)

# 筛选相关性>=0.3的模块
cormax <- cor_p_count %>% filter(abs(cor_p_count$corvalue) >= 0.3)

# 相关性>=0.3的模块的基因总数量
cormaxgenescount <- sum(cormax$number)

# 找出与下一阶梯数值相差大的那几个模块
# 增加相关性绝对值列和差值列
cor_abs_order <- cormax %>%
    arrange(-abs(corvalue)) %>%
    mutate(corabs = abs(corvalue)) %>%
    mutate(subvalue = 0)

# 计算差值
for (i in 1:(nrow(cor_abs_order) - 1)) {
    cor_abs_order$subvalue[i] <- cor_abs_order$corabs[i] - cor_abs_order$corabs[i + 1]
}

# 取梯度最大的那一批模块
select_mod0 <- cor_abs_order[1:{
    which(cor_abs_order$subvalue == max(cor_abs_order$subvalue))
}, ] %>% pull(module)

# 所挑选出来的模块的相关性：
select_modcor <- cor_abs_order[1:{
    which(cor_abs_order$subvalue == max(cor_abs_order$subvalue))
}, ] %>% pull(corvalue)
# 如果存在相关性>0.3的模块
if (nrow(cormax) > 0) {
    # 将模块按照相关性大小降序排列
    cor_abs_order <- cormax %>%
        arrange(-abs(corvalue)) %>%
        mutate(corabs = abs(corvalue)) %>%
        mutate(subvalue = 0)

    # 找出与下一阶梯数值相差大的那几个模块
    for (i in 1:(nrow(cor_abs_order) - 1)) {
        cor_abs_order$subvalue[i] <- cor_abs_order$corabs[i] - cor_abs_order$corabs[i + 1]
    }

    select_modname <- cor_abs_order[1:order(-cor_abs_order$subvalue)[1], ] %>% pull("module")

    select_modgenes <- cor_abs_order[1:order(-cor_abs_order$subvalue)[1], ] %>% pull("count")

    print(str_c('第一梯队数值的基因数 > 800，基于MM、GS选择第一梯队模块基因进行分析...'))

    # 统计模块总数：
    select_genesum <- cor_abs_order[1:order(-cor_abs_order$subvalue)[1], ] %>%
        pull(count) %>% sum()

    print(str_c("挑选的模块总数为:", length(select_modname), "    基因总数为:", select_genesum))

    # 当挑选模块的基因总数大于1000时，则进行MM&GS过滤，以下为生成MM&GS值：
    if (select_genesum > 800) {
        # merged_infor <- wgcna_net$merged_infor
        merged_infor <- read.delim(paste0(od, "res/merged_infor.txt"), check.names = F) %>% rename("Target_Genes" = "cluster2")
        pheno_module_list <- list(module = select_modname, pheno = rep("Target_Genes", length(select_modname)))

        MM_GS <- lapply(1:length(select_modname), function(i) {
            select_module <- pheno_module_list$module[i]
            select_pheno <- pheno_module_list$pheno[i]
            select_data <- merged_infor %>% filter(module %in% select_module)

            mm2 <- quantile(abs(select_data[, select_module]))[3]
            gs2 <- quantile(abs(select_data[, select_pheno]))[3]

            mm3 <- quantile(abs(select_data[, select_module]))[4]
            gs3 <- quantile(abs(select_data[, select_pheno]))[4]


            return(data.frame("MM_B" = mm2, "GS_B" = gs2, 'MM_T' = mm3, 'GS_T' = gs3))             
        }) %>%
            Filter(Negate(is.null), .) %>%
            do.call(rbind, .)
        
        MM_B <- round(quantile(MM_GS$MM_B)[3], 1) 
        GS_B <- round(quantile(MM_GS$GS_B)[3], 1) 

        MM_T <- round(quantile(MM_GS$MM_T)[3], 1) 
        GS_T <- round(quantile(MM_GS$GS_T)[3], 1) 

        # 根据生成的MM、GS值筛选hub genes, 如果用中位数得到的基因数目大于1000，则用四分位数的75%的值        
        moduleselect50 <- wgcna_select_modulegene(
            merged_infor = merged_infor, pheno_module_list = list(module = select_modname, pheno = rep("Target_Genes", length(select_modname))),
            od = paste0(od, "colormdl"), MM = MM_B, GS = GS_B
        )
        
        if(length(moduleselect50 %>% unlist()) > 800){
            moduleselect70 <- wgcna_select_modulegene(
            merged_infor = merged_infor, pheno_module_list = list(module = select_modname, pheno = rep("Target_Genes", length(select_modname))),
            od = paste0(od, "colormdl"), MM = MM_T, GS = GS_T

        )
        write.csv(moduleselect70 %>% unlist() %>% as.character(), paste0(od, "hub_gene.csv"))
        }else {
            write.csv(moduleselect50 %>% unlist() %>% as.character(), paste0(od, "hub_gene.csv"))
        }
        
    } else if(select_genesum > 50){
        print(str_c('第一梯队数值的基因数<800,选择所有相关性>0.3的第一梯队模块基因进行分析...'))        
        moduleselectall <- wgcna_select_modulegene(
            merged_infor = merged_infor, pheno_module_list = list(module = select_modname, pheno = rep("cluster2", length(select_modname))),
            od = paste0(od, "colormdl"), MM = 0, GS = 0
        )
        write.csv(moduleselectall %>% unlist() %>% as.character(), paste0(od, "hub_gene.csv"))
    }else {
       print(str_c('第一梯队数值的基因数<50,选择所有相关性>0.3的模块基因进行分析...'))
       modulegene0.3all <- wgcna_select_modulegene(
            merged_infor = merged_infor, pheno_module_list = list(module = cormax %>% pull('module'), pheno = rep("Target_Genes", length(cormax %>% pull('module')))),
            od = paste0(od, "colormdl"), MM = 0, GS = 0
        )
        write.csv(modulegene0 %>% unlist() %>% as.character(), paste0(od, "hub_gene.csv"))
    }
} else {
    print("没有相关性大于0.3的模块，进行p值显著性尝试........")
}



library(data.table)
# 根据选定参数筛选WGCNA网络模块基因（核心基因），生成Figure_{select_pheno}_{select_module}_cor.pdf
merged_infor <- fread("results/2.WGCNA/res/merged_infor.txt",data.table=F)
names(merged_infor)
p_load(WGCNA)
conflict_prefer("filter", "dplyr")
modulegene_res <- wgcna_select_modulegene(
    merged_infor = merged_infor, pheno_module_list = list(module = "turquoise", 
    pheno =rep("cluster2",1)),
    od = paste0(od, "colormdl") , MM = 0.9, GS = 0.8
)
hub_gene <- flatten_chr(modulegene_res) %>% unique()
fwrite(data.table(hub_gene=hub_gene),paste0(od,"hub_gene.csv"))


#与差异基因取交集
load("results/1.cluster/deg/deg.RData")

devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
a = list(Cluster2_deg = deg_res$DEGs, Cluster2_WGCNA = hub_gene)
p1 = ggvenn(a, c("Cluster2_deg", "Cluster2_WGCNA"), fill_color = c("blue", "red"), auto_scale = F) 
ggsave(filename = str_glue("{od}venn.pdf"))
deg_hub = intersect(deg_res$DEGs, hub_gene)

fwrite(data.table(deg_hub = deg_hub),paste0(od,"deg_hub.csv"))

load("/Pub/Users/cuiye/database/gene_map.RData")
gene_map %>% filter(gene_type == "protein_coding") %>% colnames()$gene_name
deg_hub1 = intersect(deg_hub, gene_map %>% filter(gene_type == "protein_coding") %>% pull(gene_name))


#桑基图
install.packages("networkD3")
library(networkD3)
plotdata = fread("data/dgidb_export_2022-12-26 (1).tsv")
plotdata1 = plotdata %>% select(gene, drug) %>% .[,c("drug", "gene")] %>% rename(source = drug, target = gene)

nodes <- data.frame(name=c(as.character(plotdata1$source), as.character(plotdata1$target)) %>% unique())
plotdata1$IDsource <- match(plotdata1$source, nodes$name)-1 
plotdata1$IDtarget <- match(plotdata1$target, nodes$name)-1
plotdata1$weight = 1
plotdata2 = data.frame(source = plotdata1$target %>% unique(), 
                       target = "Cluster2")
nodes[31,1] = "Cluster2"
plotdata2$IDsource <- match(plotdata2$source, nodes$name)-1 
plotdata2$IDtarget <- match(plotdata2$target, nodes$name)-1
plotdata4= table(plotdata1$target) %>% data.frame() %>% rename(source = Var1, weight = Freq)
plotdata2 = plotdata2  %>% inner_join(plotdata4)
plotdata3 = rbind(plotdata1, plotdata2)
sankeyNetwork(Links = plotdata3, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
              Value = "weight", NodeID = "name", units = 'TWh',
              fontSize = 12, nodeWidth = 30, height = 600, width = 800, sinksRight=FALSE)  


###ggalluvial 的冲击图
plotdata1 = plotdata %>% select(c("drug", "gene")) %>% mutate(cluster = "Cluster2") 
plotdata1$link <- 1
plotdata1 <- reshape::melt(plotdata1, id = 'link')

variable <- summary(plotdata1$variable)
plotdata1$flow <- rep(1:variable[1], length(variable))

#预指定颜色
mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462',
    '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
    '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5', 
    '#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b', '#BEAED4', 
    '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#4253ff', '#ff4308', '#D8D155',
    '#64495D', '#7CC767')


library(ggalluvial)

p <- ggplot(plotdata3, aes(x = variable, y = link,
    stratum = value, alluvium = flow, fill = value)) +
geom_stratum(width = 2/3) +  #冲击图中的堆叠柱形图
geom_flow(aes.flow = 'forward',width = 2/3) +  #冲击图连线绘制
scale_fill_manual(values = mycol) +  #颜色赋值
geom_text(stat = 'stratum', infer.label = TRUE, size = 4) +  #添加 lncRNA、miRNA 和 mRNA 标签
scale_x_discrete(limits = c('drug', 'gene', 'cluster')) +  #定义 lncRNA、miRNA 和 mRNA 列的展示顺序
labs(x = '', y = '') +  #去除 x 轴和 y 轴标题
theme(legend.position = 'none', panel.background = element_blank(),
    line = element_blank(), axis.text.y = element_blank(),
    axis.text.x=element_text(vjust=1,size=20))  #去除背景和图例

p

ggsave(filename = str_glue("{od}sankey diagram.pdf"), p)

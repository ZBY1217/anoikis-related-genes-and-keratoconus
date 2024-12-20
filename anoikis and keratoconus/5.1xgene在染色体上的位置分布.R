rm(list=ls())
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/UserCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
# invisible(lapply(libSources,source))
# source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)
# 导出路径
od <- str_glue("results/3.xgene_landscape/location/")
suppressWarnings(dir.create(od,recursive=TRUE))
# 导入数据
load("data/xgene1.RData")
load("data/gene_map.RData")
library(RCircos)
library(Cairo)
data(UCSC.HG38.Human.CytoBandIdeogram)
data_position <- gene_map %>%
    select(seqnames:end, gene_name) %>%
    mutate(value = 1) %>%
    filter(gene_name %in% xgene)
colnames(data_position) <- c("chromosome", "start", "end", "symbol", "value")

RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram, chr.exclude = NULL, tracks.inside = 30, tracks.outside = 0)
# chr.exclude=NULL;设置不显示的染色体，如 c(1,3)          
# tracks.inside=10;设置内部环形个数
# tracks.outside=0;设置外部环形个数  
CairoPDF(file = str_glue("{od}/Figure_xgene_location.pdf"), width = 10, height = 10, pointsize = 20, family = "Times")
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(data_position, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(data_position, name.col = 4, track.num = 2, side = "in", is.sorted = FALSE)
#第一列 染色体编号，需要与第一步导入的染色体数据一致;第二列 基因在染色体片段起始位点;第三列 基因在染色体片段结束位点;第四列 基因名
# side <- "in";指定内容在内侧的环形还是外侧的环形生成
# track.num <- 1;指定内容在第几个环形生成
# name.col <- 4;在染色体上添加基因名称， 指定内容在第几个环形生成
dev.off()

fwrite(xgene %>% data.frame(), file = "data/xgene.txt")

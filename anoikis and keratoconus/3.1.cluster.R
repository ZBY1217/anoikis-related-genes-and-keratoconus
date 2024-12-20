rm(list = ls())

## 读取包
library(pacman)
p_load(ConsensusClusterPlus, tidyverse)


#----GSE77938----

#读取数据
# load("data/clinical.RData")
load("data/expression_GSE77938.RData",verbose = T)
load("data/xgene1.RData")
#输出路径
od = "results/1.cluster"
dir.create(od, recursive = TRUE)

#表达矩阵处理
KC_sample = colnames(expression)[colnames(expression) %>% str_detect("KC")]
dat = expression[xgene, KC_sample]
#标准化
dat = log(dat+1)
##归一化操作
dat = sweep(dat, 1, apply(dat, 1, median, na.rm=T))

#聚类函数
cluster = function(dat, od, maxK = 10, file_type = "png"){
#加载包
library(pacman)
p_load(ConsensusClusterPlus, tidyverse, magrittr)

#生成聚类，距离参数
distances <- c("euclidean", "pearson", "spearman", "maximum", "binary")
clusteralgs <- c("km", "pam", "hc")
distalgs <- expand.grid(clusteralgs, distances,stringsAsFactors=FALSE)
distalgs %<>% dplyr::filter(!(Var1 == "km" & Var2 != "euclidean"))

#聚类
result = pmap(distalgs, function(Var1, Var2){
        conClust <- ConsensusClusterPlus(
            as.matrix(dat),
            maxK = maxK,
            reps = 100,
            pItem = 0.8,
            pFeature = 1,
            clusterAlg = Var1, # hc,pam,km
            distance = Var2, # pearson,spearman,euclidean,binary,maximum,canberra,minkowski
            innerLinkage = "ward.D2",
            seed = 1234,
            plot = file_type,
            title = paste0(od, "/", Var1, "_", Var2),
            writeTable = FALSE
        )
        #PAC标准筛选k值
        Kvec = 2:maxK
        x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
        PAC = rep(NA,length(Kvec))
        names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

        for(i in Kvec){
            M = conClust[[i]]$consensusMatrix
            Fn = ecdf(M[lower.tri(M)])
            PAC[i-1] = Fn(x2) - Fn(x1)
            }#end for i

        # The optimal K
        optK = Kvec[which.min(PAC)]
        conClust$best_K = optK
        conClust
})
names(result) = paste(distalgs[[1]], distalgs[[2]])
return(result)}

result = cluster(dat = dat, od = od)

##hc_maximum
cluster_hc_maximum_2 = result$`hc maximum`[[2]]$consensusClass %>% 
                       data.frame() %>% rename("cluster" = ".")
save(result, file = str_glue("{od}/result.RData"))
save(cluster_hc_maximum_2, file = str_glue("{od}/cluster.RData"))

##hc_spearman
load("results/1.cluster/result.RData",verbose = T)

cluster_hc_spearman_3 = result$`hc spearman`[[3]]$consensusClass %>%  
                       data.frame() %>% rename("cluster" = ".")

save(cluster_hc_spearman_3, file = str_glue("{od}/cluster_hc_spearman_3.RData"))


#----GSE204839----



#'
v_CCluster <- function(genelist = NULL, od = NULL, exp = NULL, clin = NULL, timecol = "time",
    statuscol = "status", seed = 123456, maxK = 5, cluster_character = "Cluster",plotFormat="pdf",
    color_fun = color_fun1, xlab = "days", input_distance = NULL, input_clusteralg = NULL) {
    if (!dir.exists(od)) {
        dir.create(od)
    }
    # 过滤生存信息不全的样本
    colnames(clin)[colnames(clin) == timecol] <- "time"
    colnames(clin)[colnames(clin) == statuscol] <- "status"
    clin <- clin %>%
        mutate(time = as.numeric(time), status = as.numeric(status)) %>%
        filter(time > 0 & status != "" & status != "NA")
    # 过滤样本
    dat_exp <- exp[which(rownames(exp) %in% genelist), which(colnames(exp) %in% clin$sample)]
    clin <- clin[match(colnames(dat_exp), clin$sample), ]
    infor <- cbind.data.frame(clin, t(dat_exp))
    suppressPackageStartupMessages(library(ConsensusClusterPlus))
    suppressPackageStartupMessages(library(survival))
    suppressPackageStartupMessages(library(survminer))
    getOptK <- function(conClust, minCls = 2, maxCls = maxK) {
        # 最佳分类数
        Kvec <- minCls:maxCls
        x1 <- 0.1
        x2 <- 0.9 # threshold defining the intermediate sub-interval
        PAC <- rep(NA, length(Kvec))
        names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
        for (i in Kvec) {
            M <- conClust[[i]][["consensusMatrix"]]
            Fn <- ecdf(M[lower.tri(M)])
            PAC[i - 1] <- Fn(x2) - Fn(x1)
        }
        optK <- Kvec[which.min(PAC)]
        return(optK)
    }
    res_list <- list()

    if (is.null(input_distance)) {
        distances <- c("euclidean", "pearson", "spearman", "maximum", "binary")
    } else {
        distances <- input_distance
    }
    if (is.null(input_clusteralg)) {
        clusteralgs <- c("km", "pam", "hc")
    } else {
        clusteralgs <- input_clusteralg
    }
    distalgs <- expand.grid(clusteralgs, distances,stringsAsFactors=FALSE)
    distalgs %<>% dplyr::filter(!(Var1 == "km" & Var2 != "euclidean"))
    clusteralg <- pull(distalgs,1)
    distance <- pull(distalgs, 2)

    for (i in seq(length(distance))) {
        j <- i
        pic_km <- list()
        kmun <- 1
        # 无监督聚类
        conClust <- ConsensusClusterPlus(
            as.matrix(dat_exp),
            maxK = maxK,
            reps = 100,
            pItem = 0.8,
            pFeature = 1,
            clusterAlg = clusteralg[j], # hc,pam,km
            distance = distance[i], # pearson,spearman,euclidean,binary,maximum,canberra,minkowski
            innerLinkage = "ward.D2",
            seed = seed,
            plot = plotFormat,
            title = paste0(od, "/", clusteralg[j], "_", distance[i]),
            writeTable = FALSE
        )
        # 计算2-maxK中最优分组的k值
        bestk <- getOptK(conClust)
        # 循环输出最优分组、2-maxK其它分组的预后曲线
        for (k in c(bestk, setdiff(2:maxK, bestk))) {
            cluster <- data.frame(
                sample = as.character(names(conClust[[k]]$consensusClass)),
                Cluster = str_c(cluster_character, conClust[[k]]$consensusClass)
            )
            dat <- data.frame(
                time = as.numeric(infor$time),
                status = as.numeric(infor$status),
                Cluster = cluster$Cluster
            ) %>% arrange(`Cluster`)
            kmplot <- survfit(Surv(time, status) ~ Cluster, data = dat)
            surv_diff <- survdiff(Surv(time, status) ~ Cluster, data = dat)
            p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
            if (length(unique(dat$Cluster)) > 3) {
                fontsize <- 2
            } else {
                fontsize <- 4
            }
            p <- ggsurvplot(kmplot,
                data = dat, conf.int = FALSE, pval = TRUE, conf.int.style = "step",
                risk.table = "absolute",
                pval.size = 5, palette = color_fun,
                legend.title = cluster_character,
                legend.labs = unique(dat$Cluster),
                fontsize = fontsize,
                risk.table.y.text = FALSE, ncensor.plot = FALSE,
                xlab = xlab
            )
            if (length(unique(dat$Cluster)) > 3) {
                p <- p + guides(color = guide_legend(nrow = 2))
            }
            pic_km[[kmun]] <- p
            kmun <- kmun + 1
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["p.val"]] <- p.val
            colnames(cluster)[2] <- cluster_character
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["cluster"]] <- cluster
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["best"]] <- paste0("k", bestk)
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["kmplot"]] <- p
        }
        pdf(file = paste0(od, "/", clusteralg[j], "_", distance[i], "/km.pdf"), width = 6, height = 6)
        print(pic_km, newPage = F)
        dev.off()
    }
    # 结果整理
    res_select <- map_dfr(names(res_list), function(distance) {
        map_dfr(names(res_list[[distance]]), function(clusteralg) {
            map_dfr(names(res_list[[distance]][[clusteralg]]), function(k) {
                tmp <- cbind.data.frame(distance = distance, clusteralg = clusteralg, k = k, p.val = res_list[[distance]][[clusteralg]][[k]][["p.val"]], bestk = res_list[[distance]][[clusteralg]][[k]][["best"]])
                return(tmp)
            })
        })
    })
    return(list(res_list = res_list, res_select = res_select))
}

source("/Pub/Users/cuiye/RCodes/UserCode/newlover/plotout.R")

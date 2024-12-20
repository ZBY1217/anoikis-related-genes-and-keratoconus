
#'
v_modelscore_km_roc <- function(signature_coef = NULL, od = NULL, exp = NULL, clin = NULL, best_cut = FALSE,
    dataset = NULL, time = c(1, 3, 5), timecol = "time", statuscol = "status",
    color_fun = color_fun1, surtime_unit = 365,
    no_roc_pheatmap = FALSE, xlab = "days") {
    if (!dir.exists(od)) dir.create(od)
    # 过滤生存信息不全的样本
    colnames(clin)[colnames(clin) == timecol] <- "time"
    colnames(clin)[colnames(clin) == statuscol] <- "status"
    clin <- clin %>%
        mutate(time = as.numeric(time), status = as.numeric(status)) %>%
        filter(time > 0 & status != "" & status != "NA")
    # 判断signature是否存在，如果不存在，设置为0
    exists_sign <- intersect(signature_coef$signature, rownames(exp))
    message(sprintf("input signature num is %s,exists signature num is %s", length(signature_coef$signature), length(exists_sign)))
    if (length(exists_sign) > 0) {
        if (length(setdiff(signature_coef$signature, exists_sign)) > 0) {
            exp <- exp %>% as.data.frame()
            for (i in setdiff(signature_coef$signature, exists_sign)) {
                exp[i, ] <- 0
            }
        }
        library(survivalROC)
        library(survival)
        library(survminer)
        infor <- exp %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(., "sample") %>%
            merge(., clin, by = "sample")
        rownames(infor) <- infor$sample
        sigdata <- as.matrix(subset(infor, select = signature_coef$signature))
        Score <- sigdata %*% signature_coef$coef
        colnames(Score) <- "Score"
        write.table(Score, file = paste0(od, "/SupplementaryTable_", dataset, "_Score.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = T)
        # 生存曲线
        dat <- cbind.data.frame(time = infor$time, status = infor$status, Score = Score)
        # 分组
        if (best_cut) {
            sur.cut <- surv_cutpoint(data = dat, time = "time", event = "status", variables = "Score")
            cut <- summary(sur.cut)$cutpoint
        } else {
            cut <- median(dat[, 3])
        }
        Group <- factor(ifelse(dat[, 3] > cut, "High", "Low"), levels = c("Low", "High")) %>% as.data.frame()
        # 样本大于一组并且每组样本数目大于5个
        if (length(unique(Group[, 1])) > 1 & sum(summary(Group[, 1]) > 5) == length(unique(Group[, 1]))) {
            colnames(Group) <- "Group"
            rownames(Group) <- rownames(Score)
            write.table(Group, file = paste0(od, "/SupplementaryTable_", dataset, "_Group.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = T)
            dat <- cbind.data.frame(time = infor$time, status = infor$status, Score = Group$Group)
            data_plot <- survfit(Surv(time, status) ~ Score, data = dat)
            data_survdiff <- survdiff(Surv(time, status) ~ Score, data = dat)
            pvalue <- 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) - 1)
            HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
            lower95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] + 1 / data_survdiff$exp[2]))
            upper95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] + 1 / data_survdiff$exp[2]))
            # 生存曲线中的标签
            if (pvalue < 0.001) {
                pval <- paste0("p < 0.001", "\n", "HR = ", round(HR, 2), "\n", "95%CI (", round(lower95, 2), "-", round(upper95, 2), ")")
            } else {
                pval <- paste0("p = ", round(pvalue, 4), "\n", "HR = ", round(HR, 2), "\n", "95%CI (", round(lower95, 2), "-", round(upper95, 2), ")")
            }
            res <- ggsurvplot(data_plot,
                data = dat, conf.int = FALSE, pval = pval, conf.int.style = "step", legend = "top", legend.title = dataset,
                risk.table = "absolute", palette = color_fun, pval.size = 5,
                risk.table.y.text = FALSE, ncensor.plot = FALSE, xlab = xlab
            )
            p1 <- ggarrange(res$plot, res$table, ncol = 1, align = "v", heights = c(0.75, 0.3))
            # plotout(p = p1, od = od, num = sprintf("_%s_group_km", dataset))

            # ROC曲线
            dat <- cbind.data.frame(time = infor$time, status = infor$status, Score = Score)
            if (!is.na(HR)) {
                if (HR < 1) {
                    dat <- dat %>% mutate(status = ifelse(status > 0, 0, 1))
                }
                library(timeROC)
                roc_model <- timeROC(
                    T = dat$time,
                    delta = dat$status,
                    marker = dat$Score,
                    cause = 1,
                    weighting = "marginal",
                    times = c(
                        time[1] * surtime_unit,
                        time[2] * surtime_unit,
                        time[3] * surtime_unit
                    ),
                    ROC = TRUE,
                    iid = FALSE
                )
                roc_plot <- data.frame(
                    year1x = roc_model$FP[, 1], year1y = roc_model$TP[, 1],
                    year2x = roc_model$FP[, 2], year2y = roc_model$TP[, 2],
                    year3x = roc_model$FP[, 3], year3y = roc_model$TP[, 3]
                )
                AUC_anno_1 <- sprintf("AUC at %s year = %s", time[1], sprintf("%.3f", roc_model$AUC[[1]]))
                AUC_anno_2 <- sprintf("AUC at %s year = %s", time[2], sprintf("%.3f", roc_model$AUC[[2]]))
                AUC_anno_3 <- sprintf("AUC at %s year = %s", time[3], sprintf("%.3f", roc_model$AUC[[3]]))

                p2 <- ggplot(data = roc_plot) +
                    geom_line(aes(x = year1x, y = year1y), size = 1.2, color = "#E31A1C") +
                    geom_line(aes(x = year2x, y = year2y), size = 1.2, color = "#377EB8") +
                    geom_line(aes(x = year3x, y = year3y), size = 1.2, color = "#007947") +
                    geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
                    egg::theme_article() +
                    theme(plot.background = element_rect(fill = "white")) +
                    annotate(geom = "line", x = c(0.5, 0.54), y = .17, colour = "#E31A1C", size = 1.2) +
                    annotate("text", x = 0.55, y = .17, size = 5, label = AUC_anno_1, color = "black", hjust = "left") +
                    annotate(geom = "line", x = c(0.5, 0.54), y = .11, colour = "#377EB8", size = 1.2) +
                    annotate("text", x = 0.55, y = .11, size = 5, label = AUC_anno_2, color = "black", hjust = "left") +
                    annotate(geom = "line", x = c(0.5, 0.54), y = .05, colour = "#007947", size = 1.2) +
                    annotate("text", x = 0.55, y = .05, size = 5, label = AUC_anno_3, color = "black", hjust = "left") +
                    labs(x = "1-Specificity", y = "Sensitivity") +
                    theme(
                        axis.text.x = element_text(face = "plain", size = 12, color = "black"),
                        axis.text.y = element_text(face = "plain", size = 12, color = "black"),
                        axis.title.x = element_text(face = "plain", size = 14, color = "black"),
                        axis.title.y = element_text(face = "plain", size = 14, color = "black")
                    )
                # plotout(p = p2, od = od, num = sprintf("_%s_group_roc", dataset))

                # 风险三连图
                dat <- cbind.data.frame(time = infor$time, status = infor$status, Group = Group, Score = Score) %>%
                    arrange(Score) %>%
                    mutate(x = 1:nrow(dat), status = ifelse(status == 0, "Alive", "Death"))
                v <- as.numeric(table(dat[, 3])["Low"]) + 0.5
                h <- cut
                #
                p3 <- ggplot(data = dat) +
                    geom_point(mapping = aes(x = x, y = Score, color = Group)) +
                    theme_bw() +
                    theme(panel.grid = element_blank(), legend.position = c(0.1, 0.85)) +
                    xlab("") +
                    scale_color_manual(values = color_fun, name = NULL) +
                    geom_hline(aes(yintercept = h), colour = "#000000", linetype = "dashed") +
                    geom_vline(aes(xintercept = v), colour = "#000000", linetype = "dashed")
                # plotout(p = p3, od = od, num = sprintf("_%s_threeplot_1", dataset), w = 8, h = 4)
                #
                p4 <- ggplot(data = dat) +
                    geom_point(mapping = aes(x = x, y = time, color = status)) +
                    theme_bw() +
                    theme(panel.grid = element_blank(), legend.position = c(0.1, 0.85)) +
                    xlab("") +
                    ylab(xlab) +
                    scale_color_manual(values = color_fun, name = NULL) +
                    geom_vline(aes(xintercept = v), colour = "#000000", linetype = "dashed")
                # plotout(p = p4, od = od, num = sprintf("_%s_threeplot_2", dataset), w = 8, h = 4)
                #
                library(pheatmap)
                library(ggplotify)
                annotation_col <- as.data.frame(Group %>% dplyr::rename(Score = Group))
                ann_colors <- list(
                    Score = c("High" = color_fun[2], "Low" = color_fun[1])
                )
                col <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
                plotdata <- sigdata[order(Score), ] %>%
                    as.matrix() %>%
                    t()
                rownames(plotdata) <- colnames(sigdata)
                p <- pheatmap::pheatmap(plotdata,
                    border = F, scale = "row", cluster_rows = F, cluster_cols = F,
                    show_colnames = F, fontsize = 8, legend = TRUE, annotation_legend = FALSE,
                    annotation_col = annotation_col, annotation_colors = ann_colors,
                    color = col
                )
                dev.off()
                p5 <- ggplotify::as.ggplot(p)
                # + theme(plot.margin = unit(c(2, 0, 2, 4), "lines"))
                # plotout(p = p5, od = od, num = sprintf("_%s_threeplot_3", dataset), w = 8, h = 4)
                # 拼图
                # plot.margin 上右下左
                g1 <- cowplot::plot_grid(p1, p2, ncol = 1, scale = c(0.95, 0.95), labels = c("AUTO"), label_size = 20)
                g2 <- cowplot::plot_grid(p3 + theme(plot.margin = unit(c(1, 5, 1, 1), "lines")),
                    p4 + theme(plot.margin = unit(c(1, 5, 1, 1), "lines")),
                    p5 + theme(plot.margin = unit(c(1, 1.4, 1, 4.5), "lines")),
                    ncol = 1, labels = c("C", "D", "E"), label_size = 20
                )
                p <- cowplot::plot_grid(g1, g2, ncol = 2, labels = NA, rel_widths = c(6, 8))
                plotout(p = p, od = od, num = sprintf("_%s_model", dataset), w = 14, h = 12)
                if (no_roc_pheatmap) {
                    g1 <- p1
                    g2 <- cowplot::plot_grid(p3 + theme(plot.margin = unit(c(1, 5, 1, 1), "lines")),
                        p4 + theme(plot.margin = unit(c(1, 5, 1, 1), "lines")),
                        ncol = 1, labels = c("B", "C"), label_size = 20
                    )
                    p <- cowplot::plot_grid(g1, g2, ncol = 2, labels = c("A", NA), scale = c(0.95, 1), label_size = 20)
                    plotout(p = p, od = od, num = sprintf("_%s_model_noroc_heat", dataset), w = 12, h = 7)
                }
                #
                stat_infor <- data.frame(
                    dataset = dataset, exists_sign_num = length(exists_sign), km_pval = pvalue, HR = HR,
                    auc_1 = roc_model$AUC[[1]], auc_2 = roc_model$AUC[[2]], auc_3 = roc_model$AUC[[3]],
                    best_cut = best_cut, cut_value = cut
                )
                return(list(Score = Score, Group = Group, infor = stat_infor))
            } else {
                return(NULL)
            }
        }
    } else {
        return(NULL)
    }
}

source("/Pub/Users/cuiye/RCodes/UserCode/newlover/plotout.R")
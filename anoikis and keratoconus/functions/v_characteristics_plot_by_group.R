
#'
v_characteristics_plot_by_group <- function(characteristics_score = NULL, Group = NULL, 
    od = NULL, color_fun = color_fun1,feature2show="all",
    type = NULL, cluster_rows = FALSE, heatplot_by_scale = TRUE) {
    if (!dir.exists(od)) dir.create(od)
    group_name <- colnames(Group)[1]
    suppressPackageStartupMessages(library(ComplexHeatmap))    
    over_sam <- intersect(characteristics_score$sample, rownames(Group))
    if (length(over_sam) == 0) {
        message("The samples of intersections is 0")
    }
    characteristics_list <- characteristics_score %>%
        dplyr::select(-sample) %>%
        colnames() %>%
        sort()
    #
    characteristics_score <- characteristics_score %>%
        dplyr::filter(sample %in% over_sam) %>%
        dplyr::select(sample, all_of(characteristics_list))
    Group <- Group[over_sam, 1] %>%
        as.data.frame() %>%
        dplyr::rename(Group = 1)
    rownames(Group) <- over_sam
    suppressPackageStartupMessages(library(ggpubr))
    infor <- Group %>%
        arrange(`Group`) %>%
        rownames_to_column(var = "sample") %>%
        merge(., characteristics_score)
    sign_stat <- sapply(characteristics_list, function(x) {
        data <- data.frame(value = infor[, x], Group = infor[, "Group"])
        if (length(table(Group)) > 2) {
            p <- kruskal.test(value ~ Group, data = data, exact=FALSE)$p.value
        } else {
            p <- wilcox.test(value ~ Group, data = data, exact=FALSE)$p.value
        }
        return(p)
    })
    sign_stat <- sign_stat %>%
        as.data.frame() %>%
        dplyr::rename(pvalue = ".") %>%
        mutate(p.signif = cut(pvalue, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("****", "***", "**", "*", " ")))
    sign_stat[is.nan(sign_stat$pvalue), "p.signif"] <- " "
    write.table(sign_stat, file = sprintf("%s/SupplementaryTable_%s_sign_stat.txt", od, type), sep = "\t", col.names = TRUE, row.names = TRUE)
    # 图1：箱线图
    w_df <- tidyr::pivot_longer(infor, names_to = "signature", cols = -c(sample, Group), values_to = "value") %>%
        arrange(`Group`) %>%
        arrange(`signature`)
    # -------------------如果feature超过50则只显示差异最显著的前20个-------------------
    n_feature <- length(characteristics_list)
    if(is.numeric(feature2show) && n_feature >= feature2show){
        feature2show <- sign_stat %>% slice_min(order_by=pvalue,n=feature2show) %>% rownames()
        characteristics_list <- intersect(feature2show,characteristics_list)
    }
    if(is.character(feature2show) && length(feature2show)>1){
        characteristics_list <- intersect(feature2show,characteristics_list)
    }
    if(length(characteristics_list)==0){
        message("no common gene")
        return()
    }
    w_df %<>% dplyr::filter(signature %in% all_of(characteristics_list))
    p <- ggboxplot(w_df, x = "signature", y = "value", fill = "Group",outlier.size = 0.75) +
        stat_compare_means(aes(group = Group), label = "p.signif", na.rm = TRUE) +
        ylab("") + theme_bw(12.5) +
        theme(
            panel.grid = element_blank(), legend.title = element_blank(),
            legend.position = "top", axis.title.x.bottom = element_blank(),
            axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
            axis.text = element_text(color = "black")
        ) + scale_fill_manual(values = color_fun)

    w <- length(characteristics_list) * 0.3 + 4.6
    h <- length(characteristics_list) * 0.03 + 4.6
    plotout(od = od, num = sprintf("_%s_group_boxplot", type), w = w, h = h, p = p)

    # 图2：独立箱线图
    boxplot_piclist <- lapply(characteristics_list, function(i) {
        # data <- w_df %>% filter(signature %in% i, value >= 0) %>% mutate(logvalue = log2(value + 1))
        data <- w_df %>%
            filter(signature %in% i) %>%
            arrange(`Group`) %>%
            arrange(`signature`)
        ylab <- ifelse(nchar(i) > 30, swr(i, nwrap = 30), i)
        pic <- ggboxplot(data = data, x = "Group", y = "value", fill = "Group") +
            # pic <- ggboxplot(data = data, x = "Group", y = "logvalue", fill = "Group") +
            stat_compare_means(na.rm = TRUE) +
            xlab("") + ylab(ylab) + theme_bw() +
            # xlab("") + ylab(str_glue("log2({i}+1)")) + theme_bw() +
            theme(
                panel.grid = element_blank(), legend.title = element_blank(), legend.position = "top"
            ) + scale_fill_manual(values = color_fun)
        if (length(unique(data$Group)) >= 4) {
            pic <- pic + guides(fill = guide_legend(nrow = 2))
        }
        return(pic)
    })
    nrow_used <- case_when(
        between(length(characteristics_list), 1, 4) ~ 1,
        between(length(characteristics_list), 5, 8) ~ 2,
        between(length(characteristics_list), 9, 12) ~ 3,
        between(length(characteristics_list), 13, 16) ~ 4,
        between(length(characteristics_list), 17, 20) ~ 5,
        length(characteristics_list) >= 21 ~ 6
    )
    w <- ceiling(length(characteristics_list) / nrow_used) * 4
    h <- nrow_used * 4
    p <- cowplot::plot_grid(plotlist = boxplot_piclist, nrow = nrow_used)
    plotout(od = od, num = sprintf("_%s_group_facet_boxplot", type), w = w, h = h, p = p)

    # 图3：热图
    characteristics_score %<>% dplyr::select(sample,all_of(characteristics_list))
    sign_stat <- sign_stat[characteristics_list,]
    # message(characteristics_list)
    col <- color_fun[1:length(unique(infor$Group))]
    names(col) <- infor$Group %>%
        unique() %>%
        as.character() %>%
        sort()
    col_anno <- HeatmapAnnotation(
        df = Group,
        col = list(
            Group = col
        ),
        # location = 1,
        annotation_label = group_name,
        show_legend = TRUE,
        annotation_name_side = "right"
    )
    row_anno <- rowAnnotation(
        "pvalue" = anno_text(
            as.matrix(sign_stat)[, "p.signif"],
            location = 0,
            gp = gpar(
                fontsize = 15
            )
        )
    )
    lgd1 <- Legend(
        # labels = c("ns", "<=0.05", "<=0.01", "<=0.001", "<=0.0001"),
        labels = c("ns", "<0.05", "<0.01", "<0.001", "<0.0001"),
        title = "Pvalue", size = 0.8,
        type = "points",
        background = "white",
        pch = c(" ", "*", "**", "***", "****"),
        nrow = 1,
        border = "black",
        grid_width = unit(1, "cm"),
        title_position = "topcenter"
    )
    lgd2 <- Legend(
        col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        at = c(-1, 0, 1),
        # direction = "horizontal",
        labels = c("<=-1", "0", ">=1")
    )
    if (heatplot_by_scale) {
        plotdata <- apply(characteristics_score %>% column_to_rownames(var = "sample"), 2, scale) %>% t() # 按列（signature）进行标准化
        colnames(plotdata) <- characteristics_score$sample
        plotdata[is.nan(plotdata)] <- 0
    } else {
        plotdata <- characteristics_score %>% column_to_rownames("sample") %>% t()
    }
    plotdata <- plotdata[, match(rownames(Group), colnames(plotdata))]
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    # rownames(plotdata) <- swr(rownames(plotdata))
    if (lapply(characteristics_list, function(x) { # 判断特征长度，如果过长，需要调整展示宽度
        return(nchar(x))
    }) %>% unlist() %>% max() >= 30) {
        p <- Heatmap(plotdata,
            cluster_columns = FALSE, cluster_rows = cluster_rows,
            top_annotation = col_anno, right_annotation = row_anno,
            show_column_dend = FALSE, show_row_dend = FALSE,
            rect_gp = gpar(col = NA),
            col = col_fun,
            column_split = Group[, 1],
            # name = "    Zscore    ",
            show_heatmap_legend = FALSE,
            row_names_side = "left",
            row_names_max_width = max_text_width(
                rownames(plotdata),
                gp = gpar(fontsize = 12)
            ),
            # heatmap_legend_param = list(direction = "horizontal"),
            show_row_names = TRUE, show_column_names = FALSE,
        ) %>%
            draw(., merge_legend = T, heatmap_legend_side = "bottom", annotation_legend_list = list(lgd1, lgd2)) %>%
            grid.grabExpr() %>%
            ggplotify::as.ggplot() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    } else {
        p <- Heatmap(plotdata,
            cluster_columns = FALSE, cluster_rows = cluster_rows,
            top_annotation = col_anno, right_annotation = row_anno,
            show_column_dend = FALSE, show_row_dend = FALSE,
            rect_gp = gpar(col = NA),
            col = col_fun,
            column_split = Group[, 1],
            # name = "    Zscore    ",
            show_heatmap_legend = FALSE,
            row_names_side = "left",
            # row_names_max_width = max_text_width(
            #     rownames(plotdata),
            #     gp = gpar(fontsize = 12)
            # ),
            # heatmap_legend_param = list(direction = "horizontal"),
            show_row_names = TRUE, show_column_names = FALSE,
        ) %>%
            draw(., merge_legend = T, heatmap_legend_side = "bottom", annotation_legend_list = list(lgd1, lgd2)) %>%
            grid.grabExpr() %>%
            ggplotify::as.ggplot() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    }
    h <- length(characteristics_list) * 0.13 + 3.7
    plotout(od = od, num = sprintf("_%s_group_heatmap", type), w = 12, h = h, p = p)
        # 图1+：violion
    p <- ggviolin(w_df, x = "signature", y = "value", fill = "Group") +
        stat_compare_means(aes(group = Group), label = "p.signif", na.rm = TRUE) +
        ylab("") + theme_bw(12.5) +
        theme(
            # panel.grid = element_blank(), 
            legend.title = element_blank(),
            legend.position = "top", axis.title.x.bottom = element_blank(),
            axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
            axis.text = element_text(color = "black")
        ) + scale_fill_manual(values = color_fun)

    w <- length(characteristics_list) * 0.3 + 4.6
    h <- length(characteristics_list) * 0.03 + 4.6
    plotout(od = od, num = sprintf("_%s_group_violin", type), w = w, h = h, p = p)

    # 图2+：独立violion
    violin_piclist <- lapply(characteristics_list, function(i) {
        # data <- w_df %>% filter(signature %in% i, value >= 0) %>% mutate(logvalue = log2(value + 1))
        data <- w_df %>%
            filter(signature %in% i) %>%
            arrange(`Group`) %>%
            arrange(`signature`)
        ylab <- ifelse(nchar(i) > 30, swr(i, nwrap = 30), i)
        pic <- ggviolin(data = data, x = "Group", y = "value", fill = "Group") +
            # pic <- ggboxplot(data = data, x = "Group", y = "logvalue", fill = "Group") +
            stat_compare_means(na.rm = TRUE) +
            xlab("") + ylab(ylab) + theme_bw() +
            # xlab("") + ylab(str_glue("log2({i}+1)")) + theme_bw() +
            theme(
                panel.grid = element_blank(), legend.title = element_blank(), legend.position = "top"
            ) + scale_fill_manual(values = color_fun)
        if (length(unique(data$Group)) >= 4) {
            pic <- pic + guides(fill = guide_legend(nrow = 2))
        }
        return(pic)
    })
    used_nrow <- case_when(
        between(length(characteristics_list), 1, 4) ~ 1,
        between(length(characteristics_list), 5, 8) ~ 2,
        between(length(characteristics_list), 9, 12) ~ 3,
        between(length(characteristics_list), 13, 16) ~ 4,
        between(length(characteristics_list), 17, 20) ~ 5,
        length(characteristics_list) >= 21 ~ 6)
    w <- ceiling(length(characteristics_list) / nrow_used) * 4
    h <- nrow_used * 4
    p <- cowplot::plot_grid(plotlist = violin_piclist, nrow = nrow_used)
    plotout(od = od, num = sprintf("_%s_group_facet_violin", type), w = w, h = h, p = p)
    return(characteristics_list)
}

source("/Pub/Users/cuiye/RCodes/UserCode/newlover/plotout.R",local = TRUE)
source("/Pub/Users/cuiye/RCodes/UserCode/newlover/swr.R",local = TRUE)
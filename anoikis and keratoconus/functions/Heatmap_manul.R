
Heatmap_manul <- function(data_input = NULL, color_used = NULL, group_infor = NULL, Colored_OtherInfor = F,
    saveplot = T, output_dir = "./", var_name = NULL, width_used = 9, height_used = 9, heatmap_name = " ",
    cluster_name = NULL, show_rownames = T, DoWilcox.test = F, rownames_fontsize = 7) {
    if (is.null(var_name)) {
        var_name <- paste0("_", paste0(sample(letters, 4), collapse = "", sep = ""))
    }

    if (!is.null(cluster_name)) {
        warning(stringr::str_glue('在function {crayon::blue("Heatmap_manul")}中，参数{crayon::bold("cluster_name")} 已弃用,并且使用{crayon::bold("group_infor")}中第一列做分组依据.'))
    }
    suppressPackageStartupMessages(library(ComplexHeatmap))
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(ggpubr))
    cluster_name <- colnames(group_infor)[1]

    sample_common <- intersect(rownames(group_infor), colnames(data_input))
    group_infor <- group_infor[sample_common, , drop = F] %>% dplyr::rename("group" = 1)
    data_input <- data_input[, sample_common]

    if (ncol(group_infor) >= 2) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")
        colnames_chr <- map_chr(clinical_names, function(x) {
            # message(x)
            # x <- "HRD"
            fisher_test_res <- group_infor %>%
                dplyr::select(group, any_of(x)) %>%
                table() %>%
                fisher.test(simulate.p.value = TRUE, B = 1e7)
            fisher_p <- fisher_test_res$p.value
            pval_label <- cut(
                x = fisher_p, breaks = c(1, .05, .01, .001, .0001, 0),
                labels = c("****", "***", "**", "*", "")
            ) %>% as.character()

            x_pval_label <- paste0(x, " ", pval_label)
        })

        colnames(group_infor)[c(1:ncol(group_infor))[-which(colnames(group_infor) %in% "group")]] <- colnames_chr
    }

    group_chara <- sort(unique(as.character(group_infor$group)))

    if (is.null(color_used)) {
        col_name <- RColorBrewer::brewer.pal(8, "Set1")[1:length(group_chara)]
        col_name <- col_name[seq_along(unique(group_infor %>% pull(group)))]
    } else {
        col_name <- color_used[seq_along(unique(group_infor %>% pull(group)))]
    }

    names(col_name) <- group_chara

    if (ncol(group_infor) >= 2) {
        col_anno <- HeatmapAnnotation(
            df = group_infor, # %>% .[, "group", drop = F]
            col = list(
                group = col_name
            ),
            annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
            annotation_name_side = "right",
            annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
        )
    } else {
        col_anno <- HeatmapAnnotation(
            df = group_infor, # %>% .[, "group", drop = F]
            col = list(
                group = col_name
            ),
            annotation_name_side = "right",
            annotation_label = cluster_name
        )
    }

    if (isTRUE(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(RColorBrewer::brewer.pal(9, "Paired"), 10)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names

        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    if (is.null(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(RColorBrewer::brewer.pal(9, "Paired"), 10)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names
        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    if (is.character(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(Colored_OtherInfor, 30)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names

        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    data_input_scaled <- t(apply(data_input, 1, function(x) scale(x, center = T, scale = T)))
    colnames(data_input_scaled) <- colnames(data_input)

    # data_input_scaled <- map_df(1:nrow(data_input), function(i) {
    #     z_scores <- (data_input[i, ] - mean(as.numeric(data_input[i, ]))) / sd(as.numeric(data_input[i, ]))
    #     return(z_scores)
    # })

    all(colnames(data_input_scaled) == group_infor %>% rownames())
    col_zscore <- circlize::colorRamp2(c(-2, 0, 2), c("#2266AC", "white", "#B2182E")) # c("#0a5aa5", "white", "firebrick3")

    # col_zscore <- circlize::colorRamp2(c(-2, 0, 2), c("#00441B", "white", "#40004B"))

    cluster_res <- Heatmap(
        matrix = data_input_scaled,
        name = heatmap_name,
        col = col_zscore,
        show_column_dend = F,
        show_column_names = F,
        show_row_names = show_rownames,
        cluster_rows = F,
        cluster_columns = F,
        show_row_dend = F,
        clustering_method_columns = "complete",
        # rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = rownames_fontsize),
        column_title = " ",
        top_annotation = col_anno,
        column_split = group_infor %>% .[, "group", drop = F],
        border = F
        # left_annotation = row_anno,
        # row_names_side = "right"
    )

    if (isTRUE(DoWilcox.test)) {
        group_wilcoxtest_p <- map_dbl(1:nrow(data_input), function(i) {
            # i <- 1
            if (length(unique(group_infor[, 1])) == 2) {
                data_used <- data.frame(value = data_input[i, ] %>% as.numeric(), group = group_infor[match(colnames(data_input), rownames(group_infor)), "group"])
                wilcox.test_res <- wilcox.test(value ~ group, data = data_used)
                pval <- wilcox.test_res$p.value

                return(pval)
            } else {
                data_used <- data.frame(value = data_input[i, ] %>% as.numeric(), group = group_infor[match(colnames(data_input), rownames(group_infor)), "group"])
                kruskal.test_res <- kruskal.test(value ~ group, data = data_used)
                pval <- kruskal.test_res$p.value

                return(pval)
            }
        })

        names(group_wilcoxtest_p) <- rownames(data_input)

        row_anno <- rowAnnotation(
            P.val = anno_text(case_when(
                between(group_wilcoxtest_p, 0.01, 0.05) ~ "*",
                between(group_wilcoxtest_p, 0.001, 0.01) ~ "**",
                between(group_wilcoxtest_p, 0.0001, 0.001) ~ "***",
                group_wilcoxtest_p < 0.0001 ~ "****",
                group_wilcoxtest_p > 0.05 ~ " "
            ),
            gp = gpar(fontsize = 8),
            location = 1,
            just = "right"
            )
        )

        cluster_res <- Heatmap(
            matrix = data_input_scaled,
            name = heatmap_name,
            col = col_zscore,
            show_column_dend = F,
            show_column_names = F,
            show_row_names = show_rownames,
            cluster_rows = F,
            cluster_columns = F,
            show_row_dend = F,
            clustering_method_columns = "complete",
            # rect_gp = gpar(col = "white", lwd = 1),
            row_names_gp = gpar(fontsize = rownames_fontsize),
            column_title = " ",
            top_annotation = col_anno,
            column_split = group_infor %>% .[, "group", drop = F],
            border = F,
            left_annotation = row_anno
            # row_names_side = "right"
        )
    }

    # gg.cluster_res <- draw(cluster_res, merge_legend = T) %>%
    #     grid.grabExpr() %>%
    #     ggplotify::as.ggplot()

    if (saveplot) {
        dir_now <- str_glue("{output_dir}")

        if (!dir.exists(dir_now)) {
            dir.create(dir_now, recursive = T)
        } else {
            message(str_c(dir_now, " is ready."))
        }

        pdf(
            file = str_glue("{dir_now}/Figure_HeatPlot_{var_name}.pdf"),
            width = width_used, height = height_used
        )
        draw(cluster_res, merge_legend = T,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()

        # ggsave2(
        #     filename = str_glue("{dir_now}/figure_heat_plot_ggsave_test{var_name}.pdf"),
        #     plot = gg.cluster_res, width = width_used, height = height_used
        # )
        # ggsave2(
        #     filename = str_glue("{dir_now}/figure_heat_plot{var_name}.tiff"),
        #     plot = gg.cluster_res, width = width_used, height = height_used, dpi = 300
        # )
        # ggsave2(
        #     filename = str_glue("{dir_now}/figure_heat_plot{var_name}_dpi72.tiff"),
        #     plot = gg.cluster_res, width = width_used, height = height_used, dpi = 72
        # )
    }
    return(cluster_res)
}
# #' @description 只对列排序，不聚类版本
# Heatmap_manul_NoSplit <- function(data_input = NULL, color_used= NULL, cluster_name= "cluster", group_infor = NULL,
#                           saveplot = T, output_dir= './', var_name= NULL, width_used = 9, height_used = 9, heatmap_name = " ") {

#     if(is.null(var_name)){
#         var_name <- paste0("_",paste0(sample(letters,4), collapse = '',sep = ""))
#     }

#     library(tidyverse)
#     library(ComplexHeatmap)
#     library(cowplot)
#     library(ggpubr)

#     sample_common <- intersect(rownames(group_infor),colnames(data_input))
#     group_infor <- group_infor[sample_common,,drop = F] %>% dplyr::select('group',everything())
#     data_input <- data_input[,sample_common]


#     group_chara <- sort(unique(as.character(group_infor$group)))

#     if (is.null(color_used)) {
#         col_name <- RColorBrewer::brewer.pal(8, "Set1")[1:length(group_chara)]
#         col_name <- col_name[seq_along(unique(group_infor %>% pull(group)))]
#     } else {
#         col_name <- color_used[seq_along(unique(group_infor %>% pull(group)))]
#     }

#     names(col_name) <- group_chara

#     if (ncol(group_infor) >= 2) {
#         col_anno <- HeatmapAnnotation(
#             df = group_infor, # %>% .[, "group", drop = F]
#             col = list(
#                 group = col_name
#             ),
#             annotation_legend_param = list(group = list(title = c(cluster_name,base::setdiff(colnames(group_infor),'group')))),
#             annotation_name_side = "right",
#             annotation_label = c(cluster_name,base::setdiff(colnames(group_infor),'group'))
#         )
#     } else {
#         col_anno <- HeatmapAnnotation(
#             df = group_infor, # %>% .[, "group", drop = F]
#             col = list(
#                 group = col_name
#             ),
#             annotation_name_side = "left",
#             annotation_label = cluster_name,
#             annotation_legend_param = list(group = list(title = cluster_name))
#         )
#     }

#     data_input_scaled <- t(apply(data_input, 1, function(x) scale(x, center = T, scale = T)))
#     colnames(data_input_scaled) <- colnames(data_input)
#     all(colnames(data_input_scaled) == group_infor %>% rownames())

#     col_zscore <- circlize::colorRamp2(c(-4, 0, 4), c("navy", "white", "firebrick3")) # c("#0a5aa5", "white", "firebrick3")

#     cluster_res <- Heatmap(
#         matrix = data_input_scaled,
#         name = heatmap_name,
#         col = col_zscore,
#         show_column_dend = F,
#         show_column_names = F,
#         show_row_names = T,
#         cluster_rows = T,
#         show_row_dend = F,
#         # rect_gp = gpar(col = "white", lwd = 1),
#         row_names_gp = gpar(fontsize = 9.5),
#         column_title = " ",
#         top_annotation = col_anno,
#         # column_order = rownames(group_infor_2),
#         column_gap = unit(0,"mm"),
#         # heatmap_legend_param = list(direction = "horizontal"),
#         column_split = group_infor %>% .[, 'group', drop = F],
#         border = F
#         # left_annotation = row_anno,
#         # row_names_side = "right"
#     )

#     gg.cluster_res <- cluster_res %>%
#         draw(., merge_legend = T,heatmap_legend_side = "right",
#                 annotation_legend_side = "right") %>%
#         grid.grabExpr() %>%
#         ggplotify::as.ggplot() +
#         theme(
#             plot.background = element_rect(fill = "white", color = "white"),
#             plot.margin = unit(c(.1, .1, .3, .1), "cm")
#         )

#     if (saveplot) {
#         dir_now <- str_glue("{output_dir}/output/heatmap/")

#         if (!dir.exists(dir_now)) {
#             dir.create(dir_now, recursive = T)
#         } else {
#             print(str_c(dir_now, " is ready."))
#         }

#         ggsave2(
#             filename = str_glue("{dir_now}/figure_heat_plot{var_name}.pdf"),
#             plot = gg.cluster_res, width = width_used, height = height_used
#         )
#         ggsave2(
#             filename = str_glue("{dir_now}/figure_heat_plot{var_name}.tiff"),
#             plot = gg.cluster_res, width = width_used, height = height_used, dpi = 300
#         )
#         ggsave2(
#             filename = str_glue("{dir_now}/figure_heat_plot{var_name}_dpi72.tiff"),
#             plot = gg.cluster_res, width = width_used, height = height_used, dpi = 72
#         )
#     }
#     return(gg.cluster_res)
# }

# 指定某几行列颜色，指定行列spilt，是否标准化，（按行列）

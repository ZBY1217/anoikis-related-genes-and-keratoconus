
#'
ClinicalFeatures_plot <- function(Data = NULL, phenotype_col = NULL, score_col = NULL, group_col = NULL,
                                    output_dir = NULL,var_name = NULL, 
                                    width_single = 3.2, height_single = 3.6,
                                    color_used = NULL, saveplot = F,
                                    savetiff = F,alpha = .8) {
    library(cowplot)
    if (is.null(var_name)) {
        var_name <- paste0("_", paste0(sample(LETTERS, 4), collapse = ""))
    }

    ignoreCase <- c(
        "NA", "NAN", "NaN", "NX", "MX", "TX", "Not Reported", "GX", "UNK", "-", "", " ",
        "unknown", "UNknown", "notreported", "gx", NA,'[Not Available]',"[Not Evaluated]"
    )

    if (!is.null(score_col)) {
        library(cowplot)
        library(tidyverse)
        library(ggpubr)

        selected_col <- c(phenotype_col, score_col) %>% unique()

        pheno_longdata <- Data %>%
            dplyr::select(any_of(selected_col)) %>%
            pivot_longer(cols = -c(all_of(score_col)), names_to = "pheno_name", values_to = "phe") %>%
            mutate(phe = as.character(phe), pheno_name = as.character(pheno_name)) %>%
            dplyr::rename(riskscore = score_col) %>% 
            as_tibble() %>%
            dplyr::filter(!phe %in% ignoreCase)

        p_list <- lapply(as.list(1:length(unique(pheno_longdata$pheno_name))), function(x) {
            set.seed(1110)
            # x <- 1
            label_y_pos <- pheno_longdata %>%
                dplyr::filter(! phe %in% ignoreCase) %>%
                dplyr::filter(pheno_name == pheno_longdata$pheno_name[x]) %>%
                dplyr::select(riskscore) %>%
                max() 

            label_y_pos <- ifelse(label_y_pos > 0 ,label_y_pos *1.08,label_y_pos * 0.9)

            label_y_pos2 <- pheno_longdata %>%
                dplyr::filter(! phe %in% ignoreCase) %>%
                dplyr::filter(pheno_name == pheno_longdata$pheno_name[x]) %>%
                dplyr::select(riskscore) %>%
                min()

            # 输出临床特征的总字符串大小，用于判断，作图时是否进行斜体标注x轴变量名
            no.pheno_chara <- pheno_longdata %>%
                dplyr::group_by(pheno_name) %>%
                summarise(n = str_count(paste0(unique(phe), collapse = "")))

            pheno_name_index <- which(no.pheno_chara$n > 20)

            if (unique(pheno_longdata$pheno_name)[x] %in% {no.pheno_chara[pheno_name_index,] %>% pull(1)}) {
                angle_used <- 15
                hjust_used <- 1
                vjust_used <- 1
            } else {
                angle_used <- 0
                hjust_used <- 0.5
                vjust_used <- 0
            }

            if (is.null(color_used)) {
                # color_used <- RColorBrewer::brewer.pal(8,'Dark2')
                color_used <- RColorBrewer::brewer.pal(8,'Set2')
            }

            p_tmp <- pheno_longdata %>%
                mutate(phe = str_replace_all(phe, pattern = " ", replacement = ""))%>%
                dplyr::filter(! phe %in% ignoreCase) %>%
                dplyr::filter(pheno_name == unique(pheno_longdata$pheno_name)[x]) %>%
                ggplot(., aes(x = factor(phe), y = riskscore)) +
                # geom_violin(aes(group = phe, fill = phe), size = .6, alpha = 0.1, width = 0.45) +
                #, fill = phe
                # geom_point(aes(colour = phe), position = position_jitter(0.12), shape = 20, alpha = 0.85) +
                geom_boxplot(aes(fill = phe),size = .5,outlier.size = .7, alpha = alpha, notch = F, width = 0.4,color = "black") +
                # aes(group = phe),
                labs(y = score_col) +
                labs(x = c(unique(pheno_longdata$pheno_name)[x])) +
                # ggsci::scale_fill_npg() +
                # ggsci::scale_color_npg() +
                scale_fill_manual(values = color_used) +
                scale_colour_manual(values = color_used) +
                egg::theme_article() +
                theme(
                    text = element_text(size = 14),
                    legend.position = "none"
                ) +
                ggpubr::stat_compare_means(aes(group = phe),label.y.npc = "top") + # label.y = label_y_pos
                # scale_y_continuous(limits = c(label_y_pos2 - 0.1, label_y_pos + 0.1)) +
                theme(
                    axis.text.x = element_text(angle = angle_used, hjust = hjust_used, vjust = vjust_used),
                    axis.text.x.bottom = element_text(color = "black")
                ) +
                theme(plot.background = element_rect(fill = "white", color = "white"))

            return(p_tmp)
        })

        names(p_list) <- unique(pheno_longdata$pheno_name)

        # 判断图片排列行数
        nrow_used <- case_when(
            between(length(p_list), 1, 4) ~ 1,
            between(length(p_list), 5, 8) ~ 2,
            between(length(p_list), 9, 12) ~ 3,
            between(length(p_list), 13, 16) ~ 4,
            between(length(p_list), 17, 20) ~ 5,
            length(p_list) > 21 ~ 6
        )

        # 根据图片排列行数 判断出图宽度
        width_used <- case_when(
            length(p_list) < 5 ~ case_when(
                length(p_list) == 1 ~ width_single * 1,
                length(p_list) == 2 ~ width_single * 2,
                length(p_list) == 3 ~ width_single * 3,
                length(p_list) == 4 ~ width_single * 4
            ),
            length(p_list) > 4 ~ width_single * 4
        )

        # 根据图片排列行数 判断出图高度
        height_used <- switch(nrow_used,
            height_single,
            height_single * 2.1,
            height_single * 3.1,
            height_single * 4.1,
            height_single * 5.1,
            height_single * 6.1
        )

        combination_plot <- plot_grid(
            plotlist = p_list,
            nrow = nrow_used,
            labels = "AUTO", label_size = 20 # label_y = 1.03
        ) +
            theme(plot.margin = unit(c(.5, .2, 0, .2), "cm"))

        if (saveplot) {
            if (!dir.exists(sprintf("%sphenotype_%s", output_dir, var_name))) {
                dir.create(sprintf("%sphenotype_%s", output_dir, var_name), recursive = T)
            }

            dir_now <- sprintf("%sphenotype_%s/", output_dir, var_name)

            # 循环导出单图
            map(1:length(p_list), function(i) {
                ggsave(
                    filename = sprintf("%sFigure_%s.pdf", dir_now, names(p_list)[i]), plot = p_list[[i]],
                    width = width_single, height = height_single
                )

                if (savetiff) {
                    ggsave(
                        filename = sprintf("%sFigure_%s.tiff", dir_now, names(p_list)[i]), plot = p_list[[i]],
                        width = width_single, height = height_single, dpi = 300
                    )
                    ggsave(
                        filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, names(p_list)[i]), plot = p_list[[i]],
                        width = width_single, height = height_single, dpi = 72
                    )
                }
         
            })

            ggsave(
                filename = sprintf("%sFigure_%s.pdf", dir_now, "phenotype_all"), plot = combination_plot,
                width = width_used, height = height_used
            )
            # ggsave(
            #     filename = sprintf("%sFigure_%s.tiff", dir_now, "phenotype_all"), plot = combination_plot,
            #     width = width_used, height = height_used, dpi = 300
            # )
            # ggsave(
            #     filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, "phenotype_all"), plot = combination_plot,
            #     width = width_used, height = height_used, dpi = 72
            # )
        }

        res_of_score_p <- list("p_list" = p_list, "combination_plot" = combination_plot, "pheno_longdata" = pheno_longdata)
    } else {
        res_of_score_p <- NULL
    }

    if (length(group_col) == 1) {
        p_list_group <- map(seq_along(unique(phenotype_col)), function(x) {

            Data_for_p <- Data %>%
                dplyr::select(any_of(c(group_col, unique(phenotype_col)[x]))) %>%
                dplyr::filter(!get0(unique(phenotype_col)[x]) %in% ignoreCase) %>%
                na.omit()

            contingency_table <- table(Data_for_p)
            fisher_test_res <- fisher.test(contingency_table, simulate.p.value = T)
            p.value <- format.pval(fisher_test_res$p.value, digits = 2L, eps = 0.001) 

            if (str_detect(string = p.value , pattern = "<")) {
                p.value  <- str_glue("Fisher.p{p.value}")
            } else {
                p.value  <- str_glue("Fisher.p={p.value}")
            }

            # 输出临床特征的总字符串大小，用于判断，作图时是否进行斜体标注x轴变量名
            no.pheno_chara <- Data[, group_col] %>%
                unique() %>%
                paste0(collapse = "") %>%
                str_count()

            if (no.pheno_chara >= 20) {
                angle_used <- 15
                hjust_used <- 1
                vjust_used <- 1
            } else {
                angle_used <- 0
                hjust_used <- 0.5
                vjust_used <- 0
            }

            # color_used默认配置
            if (is.null(color_used)) {
                # color_used <- RColorBrewer::brewer.pal(8,'Dark2')
                color_used <- RColorBrewer::brewer.pal(8,'Set2')
                # color_used <- ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10)
            }

            p_tmp <- Data_for_p %>%
                ggplot(mapping = aes_string(x = group_col, fill = unique(phenotype_col)[x])) +
                geom_bar(position = "fill", stat = "count", width = 0.7) +
                annotate("text", label = p.value, color = "black", y = 1.04, x = 1.15, size = 3) +
                labs(y = "Percent") +
                egg::theme_article(base_size = 13.5) +
                theme(legend.position = "top") + #, legend.title = element_blank()
                theme(plot.background = element_rect(fill = "white", color = "white")) +
                theme(
                    axis.text.x = element_text(angle = angle_used, hjust = hjust_used, vjust = vjust_used),
                    axis.text.x.bottom = element_text(color = "black")
                ) +
                scale_fill_manual(values = color_used)+
                theme(legend.margin = margin(.1,.1,-.3,.2,'cm'))


            if ((Data_for_p[, unique(phenotype_col)[x]] %>% unique() %>% length()) >= 3) {
                p_tmp <- p_tmp +
                    guides(fill = guide_legend(nrow = 2))
            }

            return(p_tmp)
        })

        names(p_list_group) <- unique(phenotype_col)

        # 判断图片排列行数
        nrow_used <- case_when(
            between(length(p_list_group), 1, 4) ~ 1,
            between(length(p_list_group), 5, 8) ~ 2,
            between(length(p_list_group), 9, 12) ~ 3,
            between(length(p_list_group), 13, 16) ~ 4,
            between(length(p_list_group), 17, 20) ~ 5,
            length(p_list_group) > 21 ~ 6
        )

        # 根据图片排列行数 判断出图宽度
        width_used <- case_when(
            length(p_list_group) < 5 ~ case_when(
                length(p_list_group) == 1 ~ width_single * 1,
                length(p_list_group) == 2 ~ width_single * 2,
                length(p_list_group) == 3 ~ width_single * 3,
                length(p_list_group) == 4 ~ width_single * 4
            ),
            length(p_list_group) > 4 ~ width_single * 4
        )

        # 根据图片排列行数 判断出图高度
        height_used <- switch(nrow_used,
            height_single,
            height_single * 2.1,
            height_single * 3.1,
            height_single * 4.1,
            height_single * 5.1,
            height_single * 6.1
        )

        combination_plot_group <- plot_grid(
            plotlist = p_list_group,
            nrow = nrow_used,
            labels = "AUTO", label_size = 20 # label_y = 1.03
        ) +
            theme(plot.margin = unit(c(.5, .2, 0, .2), "cm"))

        if (saveplot) {
            if (!dir.exists(sprintf("%s/phenotype_group_%s", output_dir, var_name))) {
                dir.create(sprintf("%s/phenotype_group_%s", output_dir, var_name), recursive = T)
            }

            dir_now <- sprintf("%s/phenotype_group_%s/", output_dir, var_name)

            # 循环导出单图
            map(1:length(p_list_group), function(i) {
                ggsave(
                    filename = sprintf("%sFigure_%s.pdf", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                    width = width_single, height = height_single
                )

                if (savetiff) {
                    ggsave(
                        filename = sprintf("%sFigure_%s.tiff", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                        width = width_single, height = height_single, dpi = 300
                    )
                    ggsave(
                        filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                        width = width_single, height = height_single, dpi = 72
                    )
                }


            })

            ggsave(
                filename = sprintf("%sFigure_%s.pdf", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used
            )
            ggsave(
                filename = sprintf("%sFigure_%s.tiff", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used, dpi = 300
            )
            ggsave(
                filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used, dpi = 72
            )
        }


        Clinical_in_Group <- list("p_list_group" = p_list_group, "combination_plot" = combination_plot_group)
    } else {
        Clinical_in_Group <- NULL
    }

    if (length(group_col) > 1 & length(phenotype_col) == 1) {
        p_list_group <- map(1:length(unique(group_col)), function(x) {
            Data_for_p <- Data %>%
                dplyr::select(any_of(c(phenotype_col, unique(group_col)[x]))) %>%
                dplyr::filter(!get0(unique(phenotype_col)[x]) %in% ignoreCase) %>% 
                na.omit()

            contingency_table <- table(Data_for_p)
            fisher_test_res <- fisher.test(contingency_table, simulate.p.value = T)
            p.value <- format.pval(fisher_test_res$p.value, digits = 2L, eps = 0.001)

            if (str_detect(string = p.value , pattern = "<")) {
                p.value  <- str_glue("Fisher.p{p.value}")
            } else {
                p.value  <- str_glue("Fisher.p={p.value}")
            }

            # 输出临床特征的总字符串大小，用于判断，作图时是否进行斜体标注x轴变量名
            no.pheno_chara <- Data[, unique(group_col)[x]] %>%
                unique() %>%
                paste0(collapse = "") %>%
                str_count()

            if (no.pheno_chara > 20) {
                angle_used <- 15
                hjust_used <- 1
                vjust_used <- 1
            } else {
                angle_used <- 0
                hjust_used <- 0.5
                vjust_used <- 0
            }

            # color_used默认配置
            if (is.null(color_used)) {
                color_used <- RColorBrewer::brewer.pal(8,'Set2')
                # color_used <- ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10)
            }

            p_tmp <- Data_for_p %>%
                ggplot(mapping = aes_string(x = unique(group_col)[x], fill = phenotype_col)) +
                geom_bar(position = "fill", stat = "count", width = 0.7) +
                annotate("text", label = p.value, color = "black", y = 1.04, x = 1.15, size = 3) +
                labs(y = "Percent") +
                egg::theme_article(base_size = 13.5) +
                theme(legend.position = "top") + #, legend.title = element_blank()
                theme(plot.background = element_rect(fill = "white", color = "white")) +
                theme(
                    axis.text.x = element_text(angle = angle_used, hjust = hjust_used, vjust = vjust_used),
                    axis.text.x.bottom = element_text(color = "black")
                ) +
                scale_fill_manual(values = color_used)+
                theme(legend.margin = margin(.1,.1,-.3,.2,'cm'))

            if ((Data_for_p[, unique(phenotype_col)] %>% unique() %>% length()) >= 3) {
                p_tmp <- p_tmp +
                    guides(fill = guide_legend(nrow = 2))
            }

            return(p_tmp)
        })

        names(p_list_group) <- unique(phenotype_col)

        # 判断图片排列行数
        nrow_used <- case_when(
            between(length(p_list_group), 1, 4) ~ 1,
            between(length(p_list_group), 5, 8) ~ 2,
            between(length(p_list_group), 9, 12) ~ 3,
            between(length(p_list_group), 13, 16) ~ 4,
            between(length(p_list_group), 17, 20) ~ 5,
            length(p_list_group) > 21 ~ 6
        )

        # 根据图片排列行数 判断出图宽度
        width_used <- case_when(
            length(p_list_group) < 5 ~ case_when(
                length(p_list_group) == 1 ~ width_single * 1,
                length(p_list_group) == 2 ~ width_single * 2,
                length(p_list_group) == 3 ~ width_single * 3,
                length(p_list_group) == 4 ~ width_single * 4
            ),
            length(p_list_group) > 4 ~ width_single * 4
        )

        # 根据图片排列行数 判断出图高度
        height_used <- switch(nrow_used,
            height_single,
            height_single * 2.1,
            height_single * 3.1,
            height_single * 4.1,
            height_single * 5.1,
            height_single * 6.1
        )

        combination_plot_group <- plot_grid(
            plotlist = p_list_group,
            nrow = nrow_used,
            labels = "AUTO", label_size = 20 # label_y = 1.03
        ) +
            theme(plot.margin = unit(c(.5, .2, 0, .2), "cm"))

        names(p_list_group) <- unique(group_col)

        if (saveplot) {
            if (!dir.exists(sprintf("%s/phenotype_group_%s", output_dir, var_name))) {
                dir.create(sprintf("%s/phenotype_group_%s", output_dir, var_name), recursive = T)
            }

            dir_now <- sprintf("%s/phenotype_group_%s/", output_dir, var_name)

            # 循环导出单图
            map(1:length(p_list_group), function(i) {
                ggsave(
                    filename = sprintf("%sFigure_%s.pdf", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                    width = width_single, height = height_single
                )

                if (savetiff) {
                    ggsave(
                        filename = sprintf("%sFigure_%s.tiff", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                        width = width_single, height = height_single, dpi = 300
                    )
                    ggsave(
                        filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, names(p_list_group)[i]), plot = p_list_group[[i]],
                        width = width_single, height = height_single, dpi = 72
                    )
                }
              
            })

            ggsave(
                filename = sprintf("%sFigure_%s.pdf", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used
            )
            ggsave(
                filename = sprintf("%sFigure_%s.tiff", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used, dpi = 300
            )
            ggsave(
                filename = sprintf("%sFigure_%s_dpi72.tiff", dir_now, "phenotype_group_all"), plot = combination_plot_group,
                width = width_used, height = height_used, dpi = 72
            )
        }

        res_of_group_p <- list("p_list_group" = p_list_group, "combination_plot" = combination_plot_group)
    } else {
        res_of_group_p <- NA
    }

    a <- list("Clinical_in_Score" = res_of_score_p, "Clinical_in_Group" = Clinical_in_Group, "Clinical_in_Groups" = res_of_group_p)
    return(a)
}



# test part------------

#data_12313 <- inner_join(HNSC$TCGA_OS[[2]],training[[3]]) %>% inner_join(.,tcga_subtype)

# phenotype_col12312 <- c('Age', 'Stage', 'Gender','Subtype_Selected')
# score_col12312 <- 'riskscore'

# output_dir <- "/pub/users/innertech/wyk/Project_output/tmp/"

# colnames(HNSC$TCGA_OS[[2]])

# tcga_subtype <- data.table::fread('/pub/users/innertech/wyk/Project_output/COAD_KEGG_MTOR_SIGNALING_PATHWAY/TCGA-COAD/KEGG_MTOR_SIGNALING_PATHWAY/TCGA-pancer-Subtype.20170308.tsv') %>% 
#   as_tibble() %>% 
#   dplyr::rename(sample = sampleID) %>% 
#   mutate(sample=str_c(sample,'A')) %>% 
#   left_join(training[[3]] %>% dplyr::select(sample, riskscore) ,.) %>% 
#   select(sample,Subtype_Selected)


# col_12 <- brewer.pal(5 , 'Set2')


# a <- score_in_clinical_model(Data = data_12313, phenotype_col = 'riskgroup', score_col = NULL,
#                               group_col = phenotype_col12312, 
#                              output_dir = output_dir,
#                             var_name = 'pheno_risk123', width_single = 3.6, height_single = 4.2, color_used = col_12, saveplot = T)

# 结果路径 "/pub/users/innertech/wyk/Project_output/tmp/phenotype_group_pheno_risk123"

# test part------------



#'
color_fun1 <- c("#377EB8", "#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#F781BF", "#A65628", "#FFFF33")
Box_in_One_p <- function(input, x_, y_, group_, color = color_fun1) {
    library(tidyverse)
    color <- color[seq_len(length(unique(input[[group_]])))]
    p4 <- ggplot(
        data = input,
        aes_string(x_, y_, fill = group_)
    ) +
        geom_boxplot(aes_string(col = group_), outlier.size = .2, width = .55, position = position_dodge(width = .65)) +
        scale_fill_manual(values = color) +
        scale_color_manual(values = color) +
        xlab(NULL) +
        ggpubr::theme_pubr(13) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        ggpubr::stat_compare_means(aes(group = get0(group_), label = ..p.signif..)) +
        stat_summary(fun = mean, geom = "crossbar", color = "white", width = 0.57, size = .2, position = position_dodge(width = .5))
    # ggpubr::stat_pvalue_manual(p_df,
    #     label = "p.signif", label.size = 3, bracket.size = 0, tip.length = 0.02
    # )

    return(p4)
}

#'
#' @inheritParams BoxOne_p
#' @title 药敏预测专用Segment plot
#' @param size_ 字符串，点大小对应列名
#'
SegmentPoint_p <- function(input, x_, size_, y_, color = color_fun1) {
    library(tidyverse)
    p1 <- ggplot(data = input, aes(get0(x_), reorder(get0(y_), get0(x_)))) +
        geom_segment(aes_string(xend = 0, yend = y_), linetype = 2) +
        geom_point(aes_string(size = size_), col = color[1]) +
        #   scale_size_continuous(range =c(2,8)) +
        #   scale_x_reverse(breaks = c(0, -0.3, -0.5),
        #                   expand = expansion(mult = c(0.01,.1))) + #左右留空
        ggpubr::theme_pubr() +
        labs(x = get0(size_), y = "", size = "-log10(P)") +
        theme(
            legend.position = "bottom",
            axis.ticks.y.right = element_blank(),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm")
        )
    return(p1)
}

#' @TODO 两个data.frame中，数值列计算spearman相关性
#' @title 两个data.frame中，数值列计算spearman相关性
#' @description 两个data.frame中，数值列计算spearman相关性，df1为单一数值，df2中可以有多个数值列
#' @param df1 data.frame，输入数据
#' @param value_in_df1 单一字符串，df1中用于计算相关性的列
#' @param df2 字符串，用于与df1中数值列进行相关性计算。理想状态是，除了key那一列，其余的都是数值型向量，然后与df1中value_in_df1 进行相关性计算。
#' @param key_in_df1df2 字符串，df1与df2中共有信息列
#' @param cormethod 字符串，相关性计算的方法，默认"spearman"
#' @param adjust 字符串，相关性计算的矫正P值方法，默认"BH"
#' @export
#' @return data.frame,对象
#' @author *WYK*
#'
Cor_OneAndMore <- function(df1, value_in_df1, df2, key_in_df1df2 = "sample", cormethod = "spearman", adjust = "BH") {
    require(tidyverse)
    require(psych)

    if (missing(value_in_df1)) {
        message("df1中缺少用于相关性分析的数值列，程序停止")
        return()
    }
    df1 <- df1[, c(key_in_df1df2, value_in_df1)]
    common_keys <- intersect(df2[[key_in_df1df2]], df1[[key_in_df1df2]])

    df1 <- df1[map_int(common_keys, ~ which(df1[[key_in_df1df2]] == .x)), ]
    df2 <- df2[map_int(common_keys, ~ which(df2[[key_in_df1df2]] == .x)), ]

    cor_res <- corr.test(df2 %>% select(-key_in_df1df2),
        df1 %>% select(-key_in_df1df2),
        method = cormethod, adjust = adjust
    )

    # type r p
    df3 <- inner_join(
        cor_res$r %>% as.data.frame() %>% rownames_to_column("type") %>% rename(rho = 2),
        cor_res$p %>% as.data.frame() %>% rownames_to_column("type") %>% rename(p = 2)
    )

    df2_long <- df2 %>% pivot_longer(cols = -key_in_df1df2, values_to = "value", names_to = "type")
    #             sample              type value
    #    <chr>           <chr>             <dbl>
    #  1 TCGA-06-2567-01 Camptothecin_1003 -4.21

    ComplexDF <- left_join(df1, df2_long) %>%
        left_join(df3) %>%
        as.data.frame()

    return(ComplexDF)
}


# > head(tcga_df)
#            sample status time riskgroup   riskscore HR_group
# 1 TCGA-41-2571-01      1   26       Low -0.30868822        0
# 2 TCGA-16-0846-01      1  119       Low -0.24254558        0
# 3 TCGA-06-2559-01      1  150       Low -0.14951827        0

# > IC50[1:3,1:3]
#            sample Camptothecin_1003 Vinblastine_1004
# 1 TCGA-06-2567-01         -4.212059        -5.102575
# 2 TCGA-26-5132-01         -3.802283        -4.951003
# 3 TCGA-26-5133-01         -2.793420        -6.944809

# dfx <- Cor_OneAndMore(tcga_df, "riskscore", IC50)
# top_type <- dfx %>%
#     filter(p < 0.05) %>%
#     select(type, rho) %>%
#     distinct() %>%
#     slice_max(rho, n = 10) %>%
#     pull(1)

# top_pos_durg_df <- dfx %>%
#     filter(type %in% top_type) %>%
#     rename("Drug" = "type", "logIC50" = "value") %>%
#     left_join(tcga_df %>% select(sample, riskgroup))

# a <- SegmentPoint_p(top_pos_durg_df, x_ = 'rho', size_ = 'p', y_ = 'Drug')

# b <- Box_in_One_p(input = top_pos_durg_df,
#     x_ = "Drug", y_ = "logIC50",
#     group_ = "riskgroup", color_in_p = RColorBrewer::brewer.pal(4, "Set1")[1:2]
# )

# c <- ggpubr::ggarrange(a,b)

# plotout(p = c,w = 7.4,h = 4,od = '/Pub/Users/wangyk/Project_wangyk/Codelib_YK/some_scr/DataVisualization/',name = 'aaa')

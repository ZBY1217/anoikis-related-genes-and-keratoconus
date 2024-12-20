
#'
v_uni_multi_cox <- function(od = NULL, infor = NULL, factors = c("Stage", "Age", "Gender", "Group"), dataset = NULL, coxp = 0.05,
                          timecol = "time", statuscol = "status", w = 10, h = 8) {
    if (!dir.exists(od)) dir.create(od)
    suppressPackageStartupMessages(library(forestmodel))
    suppressPackageStartupMessages(library(survival))
    suppressPackageStartupMessages(library(survminer))
    suppressPackageStartupMessages(library(tidyverse))
    # 过滤生存信息不全的样本
    colnames(infor)[colnames(infor) == timecol] <- "time"
    colnames(infor)[colnames(infor) == statuscol] <- "status"
    infor <- infor %>%
        mutate(time = as.numeric(time), status = as.numeric(status)) %>%
        filter(time > 0 & status != "" & status != "NA") %>%
        as.data.frame()
    # univariate
    common_type <- colnames(infor)[colnames(infor) %in% factors]
    # 过于字符型和因子型特征分组小于2组的
    usetype <- lapply(common_type, function(type) {
        if (length(unique(na.omit(infor[[type]]))) >= 2 && all(table(na.omit(infor[[type]])) >=2)){
            return(type)
        } else {
            return(NULL)
        }
    }) %>% unlist()
    uni_cox_model <- lapply(usetype, function(type) {
        cox <- as.formula(paste0("Surv(time, status) ~", type)) %>% coxph(data = infor)
        return(cox)
    })
    # 图片结果
    p1 <- forest_model(
        model_list = uni_cox_model,
        # panels = panels,
        # covariates = vars_for_table,
        merge_models = T,
        limits = c(-3,3),
        # return_data = T,
        recalculate_height = T,
        recalculate_width = T,
        # format_options = format_options,
        theme = theme_forest()
    )
    # 文本结果
    uni_cox_res <- map_dfr(usetype, function(type) {
        # type <- "Sex"
        cox <- as.formula(paste0("Surv(time, status) ~", type)) %>% coxph(data = infor)
        pvalue <- summary(cox)$coefficients[, "Pr(>|z|)"]
        # pvalue <- summary(cox)$wald["pvalue"]
        HR <- summary(cox)$coef[, 2]
        lower95 <- summary(cox)$conf.int[, "lower .95"]
        upper95 <- summary(cox)$conf.int[, "upper .95"]
        coef <- summary(cox)$coef[, 1]
        res <- data.frame(
            pvalue = round(pvalue, 5),
            HR = round(HR, 5),
            `Low 95%CI` = round(lower95, 5),
            `High 95%CI` = round(upper95, 5),
            `coef` = round(coef, 5)
        )
        rownames(res) <- rownames(summary(cox)$coef)
        return(res)
    })
    colnames(uni_cox_res) <- c("pvalue", "HR", "Low 95%CI", "High 95%CI", "coef")
    write.table(uni_cox_res, file = paste0(od, "/SupplementaryTable_", dataset, "_univariate_result.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    # 生成用于筛选的P值
    uni_cox_pvalue_all <- sapply(usetype, function(type) {
        cox <- as.formula(paste0("Surv(time, status) ~", type)) %>% coxph(data = infor)
        pvalue <- summary(cox)$wald["pvalue"]
        res <- data.frame(feature = type, pvalue = pvalue)
        return(res)
    }) %>%
        t() %>%
        as.data.frame()
    # 筛选用于多因素分析的特征
    sigtype <- rownames(uni_cox_pvalue_all)[which(uni_cox_pvalue_all[, "pvalue"] < coxp)]
    # multivariate
    if (length(sigtype) > 0) {
        multi_cox_model <- as.formula(paste0("Surv(time, status) ~", paste0(sigtype, collapse = "+"))) %>% coxph(data = infor)
        # 图片结果
        p2 <- forest_model(
            model = multi_cox_model,
            # panels = panels,
            # covariates = vars_for_table,
            merge_models = T,
            limits = c(-3,3),
            # return_data = T,
            recalculate_height = T,
            recalculate_width = T,
            # format_options = format_options,
            theme = theme_forest()
        )
        # 文本结果
        cox <- multi_cox_model
        pvalue <- summary(cox)$coefficients[, "Pr(>|z|)"]
        HR <- summary(cox)$coef[, 2]
        lower95 <- summary(cox)$conf.int[, "lower .95"]
        upper95 <- summary(cox)$conf.int[, "upper .95"]
        coef <- summary(cox)$coef[, 1]
        multi_cox_res <- as.data.frame(cbind(round(pvalue, 5), round(HR, 5), round(lower95, 5), round(upper95, 5), round(coef, 5)))
        colnames(multi_cox_res) <- c("pvalue", "HR", "Low 95%CI", "High 95%CI", "coef")
        write.table(multi_cox_res, file = paste0(od, "/SupplementaryTable_", dataset, "_multivariate_result.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
        p <- cowplot::plot_grid(p1, p2, ncol = 1, labels = "AUTO", label_size = 20, scale = c(0.95, 0.95))
        plotout(p = p, w = w, h = h, od = od, num = paste0("_", dataset, "_prognostic_independence"))
    } else {
        message("sig type num less than 1")
        multi_cox_res <- NULL
    }
    return(list(uni_cox_res = uni_cox_res, multi_cox_res = multi_cox_res, uni_cox_pvalue_all = uni_cox_pvalue_all))
}

source("/Pub/Users/cuiye/RCodes/UserCode/newlover/plotout.R")

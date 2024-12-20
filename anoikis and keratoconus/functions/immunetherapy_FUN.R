
#'
compute_score_in_icb <- function(exprs, model_coef = NULL) {
  exprs <- exprs[model_coef[, 1], , drop = F] %>%
    t() %>%
    as.data.frame()

  miss_gene <- model_coef[, 1][which(is.na(exprs[1, ]))]
  miss_gene_symbol <- paste0(miss_gene, collapse = ", ")

  exprs[is.na(exprs)] <- 0

  riskscore_res <- as.matrix(exprs) %*% as.matrix(model_coef[, 2, drop = F]) %>%
    as.data.frame() %>%
    dplyr::rename(riskscore = 1)

  if (length(miss_gene) == 1) {
    message(str_glue("{crayon::bgMagenta('WARNING')}\n{miss_gene_symbol} missed in exprs"))
  } else {
    message(sprintf("All ModelGenes are mapped in expression\n"))
  }

  return(riskscore_res)
}

#' @title immunetherapy_FUN 子程序，根据就不同方法计算得分
#' @description 根据输入建模方式或者自定义建模方式，来进行计算得分,也可以自定义function调用exprs、x或者其他输入参数计算得分。
#' @param model 字符串，已有的建模方法lasso,multicox,PCA1,PCA2。可以选择其中一个进行分析，也可以自定义
#' @param x 根据model中不同的方式，给到不同的东西，比如说lasso方法需要给到基因coef相关的data.frame，PCA1需要的是基因字符串向量，PCA2需要给到一个list等等
#' @param exprs 表达谱
#' @param ... 不同建模方法。传递到函数内部的function，可以调用exprs与x计算score或者直接给一个riskscore，格式有要求：样本在行名，且只有一列，第一列表头为riskscore，内容为各个样本的得分。
#' @return 只有一列的data.frame，行名为样本名称，列名为riskscore
#' @usage
#' score_res <- get_score(
#'  model = "means",
#'  x = sample(rownames(train_data$tumor_exprs)),
#'  exprs = train_data$tumor_exprs,
#'  function() {
#'    exprs[x, ] %>%
#'      t() %>%
#'      as.data.frame() %>%
#'      rowMeans() %>%
#'      as.data.frame() %>%
#'      rename(riskscore = 1)
#'  }
#' )
#' @export
#' @author *WYK*
#'
get_score <- function(model, x, exprs, ...) {
  if (!missing(...) && !{
    model %in% c("lasso", "multicox", "PCA2", "PCA1")
  }) {
    fun_defined <- function(...) {}
    body(fun_defined) <- body(match.fun(...))
    score_res <- fun_defined()
    return(score_res)
  }

  model_fun <- list(
    lasso = \(...){
      if (!is.data.frame(x)) {
        warning(str_glue("{model}要求输入x是DF，且第一列为基因名，第二列为相关系数"))
        return()
      }
      score_res <- compute_score_in_icb(exprs, model_coef = x)
      return(score_res)
    },
    multicox = \(...) {
      if (!is.data.frame(x)) {
        warning(str_glue("{model}要求输入x是DF，且第一列为基因名，第二列为相关系数"))
        return()
      }
      score_res <- compute_score_in_icb(exprs, model_coef = x)
      return(score_res)
    },
    PCA2 = \(...) {
      if (!is.list(x)) {
        warning("输入的x不是list类型，程序跳出。\nPCA2要求输入x是包含两个元素的列表，用元素1中所有基因的前两个主成分减去元素2中所有基因的前两个主成分作为模型得分。")
        return()
      }
      pacman::p_load(FactoMineR)

      x[[1]] <- intersect(x[[1]],rownames(exprs))
      x[[2]] <- intersect(x[[2]],rownames(exprs))

      if(!all(length(x[[1]])>=2 ,length(x[[2]] >=2))){
        warning("输入的基因列表中，第一个或第二个元素在表达谱中交集小于2个，程序跳出。")
        return()
      }

      pcaH <- PCA(X = exprs[x[[1]], ] %>% t(), graph = FALSE, ncp = 2)
      pcaL <- PCA(X = exprs[x[[2]], ] %>% t(), graph = FALSE, ncp = 2)
      score <- pcaH$ind$contrib[, 1] + pcaH$ind$contrib[, 2] - pcaL$ind$contrib[, 1] - pcaL$ind$contrib[, 2]
      pca_res <- data.frame(PC1PC2 = score) %>% t()

      score_res <- compute_score_in_icb(exprs = pca_res, data.frame(gene = "PC1PC2", coef = 1))
      return(score_res)
    },
    PCA1 = \(...) {
      if (!is.character(x)) {
        warning("输入的x不是character类型，程序跳出。\nPCA1要求输入x是基因集字符串向量。")
        return()
      }
      x <- intersect(x,rownames(exprs))
      pca_res <- prcomp(x = exprs[x, ], retx = F, scale = T, center = F)
      pca_res <- pca_res$rotation %>%
        as.data.frame() %>%
        mutate(PC1PC2 = PC1 + PC2) %>%
        dplyr::select(PC1PC2) %>%
        t() %>%
        as.data.frame()

      score_res <- compute_score_in_icb(exprs = pca_res, data.frame(gene = "PC1PC2", coef = 1))
      return(score_res)
    }
  )

  return(model_fun[[model]]())
}


#' @title immunetherapy_FUN 子程序，根据之前算好的score，运行KM分析
#' @description 本质上就是利用多种方法进行KM分析，调用ComplexKM_v2进行分析
#' @return list 第一个为KM曲线gg对象，第二个为计算好的df，包含score等信息
#' @author *WYK*
#'
do_KM <- function(...) {
  KM_method_in_ComplexKM <<- match.arg(KM_method,
    c("BestCut", "NormalCut", "TertileCut", "QuartileCut"),
    several.ok = TRUE
  )

  km_method_index <- list(F, F, F, F)
  names(km_method_index) <- c("BestCut", "NormalCut", "TertileCut", "QuartileCut")
  index <- map_int(KM_method_in_ComplexKM, ~ match(.x, names(km_method_index)))
  km_method_index[index] <- TRUE

  Ccomplexkm_res <- ComplexKM_v2(
    Input = clinical_for_km,
    surtime_unit = surtime_unit,
    time_col = time_col,
    status_col = status_col,
    score_col = riskscore_name,
    BestCut = km_method_index$BestCut,
    NormalCut = km_method_index$NormalCut,
    TertileCut = km_method_index$TertileCut,
    QuartileCut = km_method_index$QuartileCut,
    CutOffManual = NULL,
    GroupManul = NULL,
    legend = "top", SaveFile = F,
    od = od, w = 5, h = 5
  )

  KM_method_in_loop <- match.arg(KM_method_in_ComplexKM,
    c("BestCut_KM", "NormalCut_KM", "TertileCut_KM", "QuartileCut_KM"),
    several.ok = TRUE
  )

  km_res_filtered <- map(KM_method_in_loop, function(x) {
    km_p <- Ccomplexkm_res[[x]][[1]]
    km_group <- Ccomplexkm_res[[x]][[3]]
    return(list("KM_p" = km_p, "KM_group" = km_group))
  })

  names(km_res_filtered) <- KM_method_in_ComplexKM

  return(km_res_filtered)
}


#' @title 根据不同建模方法预测免疫治疗结果
#' @description 本function灵活考虑了多种建模方法，可以自定义riskscore计算方法
#' @param model 字符串，已有的建模方法lasso,multicox,PCA1,PCA2。可以选择其中一个进行分析，也可以自定义
#' @param x 根据model中不同的方式，给到不同的东西，比如说lasso方法需要给到基因coef相关的data.frame，PCA1需要的是基因字符串向量，PCA2需要给到一个list等等
#' @param exprs 表达谱
#' @param ... 不同建模方法。传递到函数内部的function，可以调用exprs与x计算score或者直接给一个riskscore，格式有要求：样本在行名，且只有一列，第一列表头为riskscore，内容为各个样本的得分。
#' @param clin 临床信息df
#' @param time_col 临床信息df中 时间所在列
#' @param status_col 临床信息df中，状态所在列
#' @param riskscore_name 得分名字
#' @param KM_method KM计算方法，详情见ComplexKM_v2，可以是c("BestCut_KM", "NormalCut_KM", "TertileCut_KM", "QuartileCut_KM")其中一个两个，
#' 也可以模糊匹配，比如说可以用c('N','B')代表其中的中位数分割方法，以及最优分组分割方法
#' @param color_in_plot 除KM外的配色
#' @param od 结果输出路径
#' @param w 图宽
#' @param h 图高
#' @param file_name 字符串，用于命名文件以及文件夹的字段
#' @usage
#' load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/data/forest_data.RData")
#' load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/out/2.major_part/major.RData", verbose = T)
#' load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/data/train_data.RData", verbose = T)
#' test <- immunetherapy_FUN(
#'  model = "means",
#'  x = sample(rownames(train_data$tumor_exprs)),
#'  exprs = train_data$tumor_exprs,
#'  function() {
#'    exprs[x, ] %>%
#'      t() %>%
#'      as.data.frame() %>%
#'      rowMeans() %>%
#'      as.data.frame() %>%
#'      rename(riskscore = 1)
#'  },
#'  clin = forest_data %>% mutate(Immunotherapy.Response.Group = ifelse(status ==1,'R','NoR')),
#'  time_col = "time",
#'  status_col = "status",
#'  riskscore_name = "ASS_SCORE",
#'  KM_method = c("Best",'N'),
#'  color_in_plot = NULL,
#'  savefile = T,
#'  od = "/Pub/Users/wangyk/project/Poroject/GAP220630002_LIHC_Age_gene_SC_bulk/test",
#'  w = 4,
#'  h = 4,
#'  file_name = "12344"
#' )
#'
#' @export
#' @return NULL
#' @author *WYK*
#'
immunetherapy_FUN <- function(model, x, exprs, ...,
                              clin, time_col = "time", status_col = "status", riskscore_name = "riskscore",
                              KM_method = c("Best_Cut"), color_in_plot = NULL, savefile = T, surtime_unit = 365,
                              od, w, h, file_name = "cohort") {
  # 调用脚本：
  source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/some_scr/NewLover/ComplexKM.R")
  source("src/functions/ClinicalFeatures_plot.R")

  pacman::p_load(rlang, magrittr, tidyverse)

  score_res <- get_score(model = model, x = x, exprs = exprs, ... = ...) %>%
    dplyr::rename(!!riskscore_name := "riskscore")
    print(colnames(clin))
    print(colnames(score_res))
  clinical_for_all_patients <- clin %>% full_join(rownames_to_column(score_res, "sample"))
  clinical_for_km <- clin %>% inner_join(rownames_to_column(score_res, "sample"))

  if (sum(is_in(c(time_col, status_col), colnames(clin))) == 2) {
    km_fun_inner <- \(...){}
    body(km_fun_inner) <- body(do_KM)
    km_res <- km_fun_inner()
  } else {
    message(crayon::blue("临床信息有缺失，KM分析跳过"))
  }

  if ("Immunotherapy.Response.Group" %>% is_in(colnames(clin))) {
    clinical_for_all_patients_ <- clinical_for_all_patients %>% 
      dplyr::rename("Response " = "Immunotherapy.Response.Group")
    box_in_all_patients <- ClinicalFeatures_plot(
      Data = clinical_for_all_patients_, phenotype_col = "Response ", score_col = riskscore_name,
      color_used = color_in_plot, saveplot = F, alpha = .9
    )

    score_in_all_patient <- box_in_all_patients$Clinical_in_Score[["p_list"]][["Response "]]

    if (isTRUE(savefile)) {
      od_now <- str_glue("{od}/{file_name}/")
      plotout(
        od = od_now,
        name = str_glue("Score_in_all_patients_{file_name}"), w = w * .88, h = h,
        p = score_in_all_patient
      )
    }

    if (!exists("km_res", envir = current_env())) {
      return()
    }

    os_related_plots <- map(KM_method_in_ComplexKM, function(method) {
      df <- km_res[[method]][["KM_group"]]

      if ("Response" %in% colnames(df)) {
        df %<>% rename("Response_raw" = "Response")
      }

      df %<>% rename("Response" = "Immunotherapy.Response.Group") %>% rename(!!riskscore_name := "Score")

      box_in_os_patients <- ClinicalFeatures_plot(
        Data = df, phenotype_col = "Response", score_col = riskscore_name,
        group_col = "Group",
        color_used = color_in_plot, saveplot = F, alpha = .9
      )
      score_in_os_patients_response <- box_in_os_patients$Clinical_in_Score[["p_list"]][["Response"]]
      percentage_stack_in_os_patients <- box_in_os_patients[["Clinical_in_Group"]][["p_list_group"]][["Response"]]

      res <- list(
        "score_in_os_patients" = score_in_os_patients_response,
        "percentage_stack_in_os_patients" = percentage_stack_in_os_patients
      )

      return(res)
    })
    names(os_related_plots) <- KM_method_in_ComplexKM

    if (isTRUE(savefile)) {
      walk(KM_method_in_ComplexKM, function(method) {
        od_now <- str_glue("{od}/{file_name}/{method}/")

        plotout(
          od = od_now,
          name = str_glue("KM_{file_name}"), w = w, h = h,
          p = km_res[[method]][["KM_p"]]
        )

        plotout(
          od = od_now,
          name = str_glue("Score_in_response_{file_name}"), w = w * .88, h = h,
          p = os_related_plots[[method]][["score_in_os_patients"]]
        )

        plotout(
          od = od_now,
          name = str_glue("Percentage_Stack_{file_name}"), w = w * .88, h = h,
          p = os_related_plots[[method]][["percentage_stack_in_os_patients"]]
        )

        p_ <- cowplot::plot_grid(
          km_res[[method]][["KM_p"]],
          os_related_plots[[method]][["score_in_os_patients"]],
          os_related_plots[[method]][["percentage_stack_in_os_patients"]],
          nrow = 1, labels = "AUTO"
        )

        plotout(
          od = od_now,
          name = str_glue("{file_name}"), w = w * 2.8, h = h,
          p = p_
        )



        write_tsv(
          x = km_res[[method]][["KM_group"]],
          file = str_glue("{od_now}/Table_Clinical_Score_{file_name}.tsv")
        )
      })
    }
  } else {
    message(crayon::blue("该队列反应信息有缺失，Score是否有显著分布差异以及百分比堆积图等分析跳过"))
  }

  if (!exists("km_res", envir = current_env())) {
    return()
  }

  if (isTRUE(savefile)) {
    walk(KM_method_in_ComplexKM, function(method) {
      od_now <- str_glue("{od}/{file_name}/{method}/")

      plotout(
        od = od_now,
        name = str_glue("KM_{file_name}"), w = w, h = h,
        p = km_res[[method]][["KM_p"]]
      )
      write_tsv(
        x = km_res[[method]][["KM_group"]],
        file = str_glue("{od_now}/Table_Clinical_Score_{file_name}.tsv")
      )
    })
  }

  return(km_res[[length(km_res)]])

}

# # 测试部分------
# load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/data/forest_data.RData")
# load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/out/2.major_part/major.RData", verbose = T)
# load("/Pub/Users/wangyk/project/Poroject/GDT220812001_GBM_immuneCellDeath_WGCNA/data/train_data.RData", verbose = T)
# library(YK.pkg)
# test <- immunetherapy_FUN(
#   model = "means",
#   x = sample(rownames(train_data$tumor_exprs)),
#   exprs = train_data$tumor_exprs,
#   function() {
#     exprs[x, ] %>%
#       t() %>%
#       as.data.frame() %>%
#       rowMeans() %>%
#       as.data.frame() %>%
#       rename(riskscore = 1)
#   },
#   clin = forest_data %>% mutate(Immunotherapy.Response.Group = ifelse(status == 1, "R", "NoR")),
#   time_col = "time",
#   status_col = "status",
#   riskscore_name = "ASS_SCORE",
#   KM_method = c("Best", "N"),
#   color_in_plot = NULL,
#   savefile = T,
#   od = "/Pub/Users/wangyk/project/Poroject/GAP220630002_LIHC_Age_gene_SC_bulk/test",
#   w = 4,
#   h = 4,
#   file_name = "12344"
# )

# -----

parameter_list <- list(
  # >>> 数据整理、差异基因、SNV、CNV展示需要的参数 >>> begin
  cancer = c("HNSC", "ORAL"),
  version = "v2",
  clinical_feature_list = c("Age", "Stage", "Sex"), # 用于比较基因表达差异和独立预后验证的临床分组信息
  n_exprs = 30, # 可以写 数字或"all",表示画图的基因
  n_gene_cell = 30, # 可以写 数字或"all", 表示画图的基因与免疫细胞相关性热图的基因数目
  n_snv = 30, # 数字, 表示画图的基因数目
  n_cnv = 30, # 数字, 表示画图的基因数目

  # *** 以下参数不需要修改 ***
  train_cohort = "TCGA",
  train_survival_outcome = "OS",
  valid_survival_outcome_list = c("OS"), # 验证集OS;PDF;DSS等生存
  drug_list = NULL, # 字符向量，药物名，默认NULL；如果指定的话，则只对指定药物计算和绘图
  # <<< before train <<< end

  # >>> 测试所的参数,一般只需要需改前两个参数 >>> begin
  cluster_row_number = c(2,3), # 聚类得到的filter_xgene_cc_res.txt结果中对应的cluster_row_number列
  model = "lasso", # 可以选择PCA1；PCA2；lasso；multicox

  # *** 以下参数不需要修改 ***
  seed = 141592653, # 随机数种子，一般不需要改
  maxK = 5, # 最大聚类数
  input_distance = c("euclidean", "pearson", "spearman"), # 聚类的距离计算方式，用于循环
  input_clusterAlg = c("hc", "km", "pam"), # 聚类方法，用于循环
  xgene_select_k = c("k2","k3"), # 选择的聚类数，用于循环
  gene_cox_bygroup = TRUE, # cox回归时，基因表达量是否要分高低两组
  logfc_list = c(.585, 1), # limma差异分析时，log2FC的取值，用于循环
  degp_list = c(.05, .01), # limma差异分析时，p的取值，用于循环
  coxp_list = c(.05, .01), # 单因素cox回归时，p的取值，用于循环
  roc_time = c(1, 3, 5), # 预后模型中的ROC时间1，3，5年
  cohort_list = NULL, # 字符向量，验证集队列向量，默认NULL；如果指定的话，则只对指定队列进行验证
  # <<<  parameters for filter <<< end

  # >>>  测试主干得到的最佳参数 >>> begin
  select_k = "k3", # 最后选择的参数，选择的聚类数
  select_clusterAlg = "pam", # 最后选择的参数，聚类方法
  select_distance = "pearson", # 最后选择的参数，聚类的距离计算方式
  log2fc = 1, # 最后选择的参数，limma差异分析时，log2FC的取值
  deg_p = 0.05, # 最后选择的参数，limma差异分析时，p的取值
  coxp = .05 # 最后选择的参数，单因素cox回归时，p的取值
  # <<<  final parameters <<< end
)

setwd("/Pub/Users/renxz/project/F221009004")
rm(list = ls())
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))


#数据处理
load("data/xgene1.RData",verbose = T)
load("results/1.cluster/cluster.RData",verbose = T)
# cluster_hc_maximum_2

cluster_hc_maximum_2 %<>% mutate(cluster = ifelse(cluster == 1, "C1", "C2")) %>%
    rownames_to_column("sample") %>%
    mutate(sample = str_remove(sample,"_"))

df <- readxl::read_excel("shouhou/ejhg20174x16.xlsx",sheet  = 1)
df %<>% rename(sample = 1)  %>% na.omit()
dim(cluster_hc_maximum_2)


df %<>%
    filter(sample %in% (cluster_hc_maximum_2$sample %>% str_remove("_"))) %>%
    as.data.frame() %>%
    rename(Age = 4, Eye_rubbing = `Eye rubbing`, Smoking = "Smok ing")

df %<>% inner_join(cluster_hc_maximum_2)
df %<>% mutate(Age = ifelse(Age > 40, ">40", "<=40"))

source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/some_scr/Clinical_features.R")


a <- score_in_clinical_model(
    Data = df, phenotype_col = c("Gender", "Age"), score_col = NULL, group_col = "cluster",
    output_dir = "shouhou/", var_name = "pheno_in_cluster",
    width_single = 3.4, height_single = 3.8,
    color_used = NULL, saveplot = T,
    savetiff = F
)


# GSE151631 数据处理
files <- list.files("shouhou/GSE151631/expr",full.names = T)

files_name <- list.files("shouhou/GSE151631/expr") %>%
    str_remove(".genes.fpkm_tracking.gz") %>%
    str_split("_") %>%
    map(~ .x[1]) %>%
    unlist()

l <- map(files, ~ data.table::fread(.x) %>% select(gene_id, FPKM))

l2 <- map2(l,files_name, ~ rename(.x, {{.y}}:= 2))
l3 <- map(l2, ~ {
    # .x = l2[[1]]
    gene <- .x[[1]]
    sample_name <- colnames(.x)[2]

    m <- .x[[2]] %>% as.matrix()
    colnames(m) <- sample_name
    rownames(m) <- gene
    m
})
gene <- rownames(l3[[1]])

l4 <- map_dfc(l3, ~ .x[gene,,drop = F])
l4 %<>% mutate(gene = gene) %>% select(gene,everything())

source("/Pub/Users/cuiye/RCodes/UserCode/newlover/exp_martix_aggregate.R")
expr <- exp_martix_aggregate(input_exp = as.data.frame(l4), rownames_type = c("gene"),fun = "median")
expr <- log2(expr + 1)

df <- read.delim("shouhou/GSE151631/clinical.txt") %>%
    rename(sample = 1) %>%
    mutate(group = ifelse(Disease == "healthy control", "control", "case")) %>%
    select(sample, group, Age, Sex,Atopy.allergy ) %>% 
    mutate(Age = ifelse(Age < 40,'<40',">=40"))


GSE151631 = list(data_exprs = expr,data_clinical = df)
saveRDS(GSE151631,'shouhou/GSE151631/GSE151631.rds')

GSE151631$data_exprs[xgene,]

source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/clucter_degs_major/cluster_degs_major_function.R")

fake_clin <- GSE151631$data_clinical %>% mutate(
    time = sample(200:1500,
        size = nrow(GSE151631$data_clinical)
    ),
    status = sample(0:1, size = nrow(GSE151631$data_clinical),replace = T)
)



res <- CC_cluster(
    exp = GSE151631$data_exprs[xgene, fake_clin %>% filter(group == 'case') %>% pull(1)],
    clinical = fake_clin %>% filter(group == 'case'), seed = 1110, k.max = 5, 
    output_dir = "/Pub/Users/renxz/project/F221009004/shouhou/"
)

cc_df <- data.table::fread("shouhou/output/CC_cluster/CC_km_euclidean/untitled_consensus_cluster/untitled_consensus_cluster.k=2.consensusClass.csv") %>%
    rename(sample = 1, cluster = 2) %>%
    mutate(cluster = ifelse(cluster == 1, "C1", "C2"))

df <- inner_join(cc_df,GSE151631$data_clinical)
source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/some_scr/Clinical_features.R")

a <- score_in_clinical_model(
    Data = df, 
    phenotype_col = c("Sex", "Age",'Atopy.allergy'),
    score_col = NULL,
    group_col = "cluster",
    output_dir = "shouhou/",
    var_name = "pheno_in_cluster_GSE151631",
    width_single = 3.4,
    height_single = 3.8,
    color_used = NULL,
    saveplot = T,
    savetiff = F
)



# GSE204839--------
expr <- read_delim("data/GSE204839GSE204839_GPL21185/raw_expression.tsv")
expr %<>% as.data.frame() %>% column_to_rownames("Gene_Symbol")

clin <- read_tsv("data/GSE204839GSE204839_GPL21185/raw_clinical.tsv")
clin %<>%
    select(geo_accession, characteristics_ch1.1, characteristics_ch1.2) %>%
    rename(sample = 1, group = 3, Sex = 2) %>%
    mutate(
        group = ifelse(group == "disease: healthy control", "control", "case"),
        Sex = ifelse(Sex == "gender: F", "F", "M")
    )
GSE204839 <- list(data_expr = expr, data_clinical = clin)
saveRDS(GSE204839, "shouhou/GSE204839.rds")


fake_clin <- GSE204839$data_clinical %>% mutate(
    time = sample(200:1500,
        size = nrow(GSE204839$data_clinical)
    ),
    status = sample(0:1, size = nrow(GSE204839$data_clinical),replace = T)
) %>% as.data.frame()

source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/clucter_degs_major/cluster_degs_major_function.R")
res <- CC_cluster(
    exp = GSE204839$data_expr[xgene, fake_clin %>% filter(group == 'case') %>% pull(1)],
    clinical = fake_clin %>% filter(group == 'case'), seed = 1110, k.max = 5, 
    output_dir = "/Pub/Users/renxz/project/F221009004/shouhou/GSE204839/"
)
cc_df <- data.table::fread("shouhou/GSE204839/output/CC_cluster/CC_hc_spearman/untitled_consensus_cluster/untitled_consensus_cluster.k=2.consensusClass.csv") %>%
    rename(sample = 1, cluster = 2) %>%
    mutate(cluster = ifelse(cluster == 1, "C1", "C2"))

df <- inner_join(cc_df,GSE204839$data_clinical)
source("/Pub/Users/wangyk/Project_wangyk/Codelib_YK/some_scr/Clinical_features.R")
a <- score_in_clinical_model(
    Data = df, 
    phenotype_col = c("Sex"),
    score_col = NULL,
    group_col = "cluster",
    output_dir = "shouhou/",
    var_name = "pheno_in_cluster_GSE204839",
    width_single = 3.4,
    height_single = 3.8,
    color_used = NULL,
    saveplot = T,
    savetiff = F
)



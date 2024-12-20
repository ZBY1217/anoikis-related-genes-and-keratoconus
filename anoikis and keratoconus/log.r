GSE[2]$platform_id[1]
View(GSE[2])

#表达谱
dim(expression)
colnames(expression)
expression$"Gene biotype" %>% table()
head(expression)

#cluster
colnames(expression)
View(clinical)
colnames(clinical)
xgene
dim(expression)
head(dat)
dim(dat)

colnames(expression)
range(dat)
range(expression)

##pca&&deg
range(expression)
dim(expression)
rownames(expression)

##富集
View(geneSets[[2]])
geneSets[[2]] @ geneIds
geneSets[[1]]@setName
length(geneSets)
head(aa)
View(gsva_matrix)
dim(dat)
colnames(gsva_matrix)
View(kegg)
head(path)


#免疫浸润
rownames(expression) %>% 
str_detect("^HLA") %>% table()
expression %>% colnames()
chemokine %in% rownames(expression)
library(data.table)
#deg
dat %>% head()

#wgcna
cluster_hc_maximum_2
library(tidyverse)
colnames(dat)
sample_data %>% head()
dat %>% head()
range(datExpr)
datExpr %>% head()
od


#相关性
dim(expression)
colnames(expression)
colnames(expression_disease)
cluster_hc_maximum_2 %>% dim()
dim(expression_xgene)
library(RColorBrewer)
display.brewer.all(type = "div")
range(cor(t(expression_disease_xgene)))
range(expression_disease_xgene)



#xgene不同分组间表达情况
expression %>% dim()
range(expression)


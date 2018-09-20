library(tidyverse)
study_gene <- readRDS("./data/processed/study_geneMat.rds")
study_drug <- readRDS("./data/processed/study_drugMat.rds")
study_pathway <- readRDS("./data/processed/study_pathwayMat.rds")
genage_genes <- readRDS("./data/processed/lit_geneMat.rds")
drugage_mat <- readRDS("./data/processed/lit_drugMat.rds")
genage_pathway <- readRDS("./data/processed/lit_pathwayMat.rds")

drugs_by_studies <- names(which(colSums(study_drug==1)>=1))
drugageDrugs <- names(which(drugage_mat[3,]==1))

# How many drugs in DrugAge
length(unique(drugageDrugs))
# [1] 346
# How many are discovered by at least one study
sum(drugageDrugs %in% drugs_by_studies)
# [1] 41
mean(drugageDrugs %in% drugs_by_studies)
# [1] 0.1184971 

# How many drugs per each study 
as.data.frame(rowSums(study_drug))
# Aliper2016                       7
# Barardo2017                     11
# Calvert2016                      7
# Donertas2018                    17
# Fernandes2016                   14
# Fuentealba2018                  10
# Liu2016                         54
# Mofidifar2018                    4
# Snell2016                       21
# Snell2018                       23
# Yang2018                         2
# Ziehm2017                       11
#on average
mean(rowSums(study_drug))
# [1] 15.08333

# How many drugs in total
length(unique(drugs_by_studies))
# 163

# How many are discovered by only one study
sum(colSums(study_drug==1)==1)
# [1] 149
mean(colSums(study_drug==1)==1)
# [1] 0.9141104

as.data.frame(table(colSums(study_drug==1))) %>%
    rename(numStudy = Var1)
# numStudy Freq
#        1  149
#        2   10
#        3    4

# How many drugs not in drugAge
sum(!drugs_by_studies %in% drugageDrugs)
# [1] 122

drugTargetsStat <- read_tsv('data/processed/study_drug_gene_literature_masterList.tsv') %>%
    filter(is.na(Drug_in_DrugAge)) %>%
    group_by(Drug, DrugName) %>%
    summarise(nGene= length(unique(Gene)), 
              Gene_in_GenAge_Human=any(Gene_in_GenAge_Human,na.rm = T),
              Gene_in_GenAge_Model=any(Gene_in_GenAge_Model,na.rm = T)) %>%
    mutate(inGenAge = Gene_in_GenAge_Human | Gene_in_GenAge_Model)

# how many novel drugs    
mean(drugTargetsStat$inGenAge)
# 0.6639344
sum(drugTargetsStat$inGenAge)
# [1] 81

genagehuman = unique(names(which(genage_genes[1,]==1)))
genagemodel = unique(names(which(genage_genes[2,]==1)))
studygenes = unique(names(which(colSums(study_gene==1)>=1)))

# How many genage human genes targeted by the drugs found in the study
length(genagehuman)
# [1] 307
sum(genagehuman %in% studygenes)
# [1] 103
mean(genagehuman %in% studygenes)
# [1] 0.3355049


# How many genage model organism genes targeted by the drugs found in the study
length(genagemodel)
# [1] 902
sum(genagemodel %in% studygenes)
# [1] 94
mean(genagemodel %in% studygenes)
# [1] 0.1042129

sort(colSums(study_gene[,intersect(names(which(colSums(study_gene)>1)),genagehuman)]),dec=T)
# DDIT3    ERBB2      BAX     BCL2     BDNF     EGFR      MYC   PIK3CA      RB1 
# 8        8        7        7        7        7        7        7        7 
# TP53 
# 7

sort(colSums(study_gene),dec=T)[1:10]
# ABCB1 BIRC5  KRAS DDIT3 ERBB2   APC   BAX  BCL2  BDNF  BRAF 
# 10     9     9     8     8     7     7     7     7     7 

read_tsv('data/processed/study_drug_gene_literature_masterList.tsv') %>%
    group_by(Drug, DrugName) %>%
    summarise(nGene= length(unique(Gene)), 
              Gene_in_GenAge_Human=all(Gene_in_GenAge_Human),
              Gene_in_GenAge_Model=all(Gene_in_GenAge_Model)) %>%
    mutate(inGenAge = Gene_in_GenAge_Human | Gene_in_GenAge_Model) %>%
    summary()

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                        ChEMBLID) %>% unique()

DA_drugs_incSyn_CHEMBLlist <- read_tsv("./data/processed/drugAgeCHEMBL.tsv")

drugageInterations <- DA_drugs_incSyn_CHEMBLlist %>%
    filter(!is.na(ChEMBLID)) %>%
    left_join(interactions) %>%
    unique()

drugageInterations2 <- drugageInterations %>%
    group_by(ChEMBLID) %>%
    filter(!is.na(gene_name)) %>%
    summarise(any(gene_name %in% studygenes))

# drugage targeting genes discovered by the studies
mean(drugageInterations2$`any(gene_name %in% studygenes)`)
# [1] 0.8026316

# pathways not targeted by the studies
sum(colSums(study_pathway!=0)==0)
# 25
# pathways targeted by the studies
sum(colSums(study_pathway!=0)>0)
# 294
mean(colSums(study_pathway!=0)>0)
# [1] 0.9216301


# pathways without genage_human genes
sum(genage_pathway[1,]==0)
# 84
# with
sum(genage_pathway[1,]>0)
# 235
mean(genage_pathway[1,]>0)
# [1] 0.7366771

# pathways without genage_model genes
sum(genage_pathway[2,]==0)
# 54
# with
sum(genage_pathway[2,]>0)
# 265
mean(genage_pathway[2,]>0)
# [1] 0.830721

# pathways without genes targeted by drugage drugs
sum(genage_pathway[3,]==0)
# 37
# with
sum(genage_pathway[3,]>0)
# 282
mean(genage_pathway[3,]>0)
# [1] 0.8840125

names(which(colMeans(study_pathway!=0)==1))
# [1] "hsa04932" "hsa05160" "hsa05200" "hsa05202"

study_pathway[,names(which(colMeans(study_pathway!=0)==1))]

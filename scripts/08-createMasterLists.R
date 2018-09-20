library(tidyverse)

study_pathway <- readRDS("./data/processed/study_pathwayMat.rds")
genage_genes <- readRDS("./data/processed/lit_geneMat.rds")
drugage_mat <- readRDS("./data/processed/lit_drugMat.rds")
genage_pathway <- readRDS("./data/processed/lit_pathwayMat.rds")

genage_genes <- reshape2::melt(genage_genes) %>% 
    filter(value == 1) %>% 
    select(-value) %>% 
    rename(Gene = Var2) %>% 
    mutate(value = T) %>% 
    spread(Var1, value, fill = NA) %>% 
    rename(Gene_in_GenAge_Human = Human, Gene_in_GenAge_Model = Model, 
           Gene_targeted_by_DrugAge = DrugAge) %>% 
    unique()

genage_drugs <- reshape2::melt(drugage_mat) %>% 
    filter(value == 1) %>% 
    select(-value) %>% 
    rename(Drug = Var2) %>% 
    mutate(value = T) %>% 
    spread(Var1, value, fill = NA) %>% 
    rename(Drug_targets_GenAge_Human = Human, Drug_targets_GenAge_Model = Model, 
           Drug_in_DrugAge = DrugAge) %>% unique()

combinedList <- read_tsv("./data/processed/combinedList.tsv") %>% 
    select(-interaction_claim_source) %>% 
    rename(Drug = ChEMBLID, Study = study, 
           Gene = gene_name, DrugName = drug_name) %>% 
    unique() %>% 
    select(Study, Drug, DrugName, Gene)

combinedList <- left_join(combinedList, genage_genes) %>% 
    left_join(genage_drugs) %>% 
    unique()

write_tsv(combinedList, "./data/processed/study_drug_gene_literature_masterList.tsv")


pathNames <- sapply(unique(union(colnames(study_pathway), colnames(genage_pathway))), 
                    function(kg) {
                        strsplit(KEGGREST::keggGet(kg)[[1]]$NAME, " - ")[[1]][1]
                    })

study_pathway <- reshape2::melt(study_pathway) %>% 
    rename(Study = Var1, pathID = Var2, percentageTargeted = value) %>% 
    mutate(Pathway = pathNames[pathID]) %>% 
    select(Study, pathID, Pathway, percentageTargeted)

genage_pathway <- reshape2::melt(genage_pathway) %>% 
    rename(pathID = Var2) %>% 
    spread(Var1, value) %>% 
    rename(percentage_in_Human_GenAge = Human, 
           percentage_in_Model_GenAge = Model, 
           percentage_targeted_by_DrugAge = DrugAge)

combinedPathway <- left_join(study_pathway, genage_pathway) %>% 
    unique()

write_tsv(combinedPathway, "./data/processed/study_pathway_literature_masterList.tsv")

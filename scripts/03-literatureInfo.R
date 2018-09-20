# libraries
library(tidyverse)

# functions
name2CID <- function(nm) {
    library(RCurl)
    library(jsonlite)
    nm <- URLencode(nm, reserved = T)
    name2cid <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", 
                             nm, "/cids/JSON", sep = ""))
    name2cid <- fromJSON(name2cid)
    return(name2cid$IdentifierList$CID)
}

convertIDs <- function(drugList) {
    cids <- sapply(drugList, name2CID)
    cids <- cids[!sapply(cids, is.null)]
    cids <- reshape2::melt(cids) %>% rename(CID = value, 
                                            drugName = L1)
    noCID <- setdiff(drugList, cids$drugName)
    return(list(cids, noCID))
}

cid2chembl <- function(intab, out_fx) {
    pubchem2chembl <- read_tsv("./data/raw/src1src22.txt") %>% 
        setNames(., c("ChEMBLID", "CID")) %>% mutate(CID = as.character(CID))
    fintable <- intab %>% mutate(CID = as.character(CID)) %>% 
        left_join(pubchem2chembl) %>% unique()
    write_tsv(fintable, out_fx)
    return(fintable)
}

# convert drugage drugs to 1)pubchem 2)chembl ids
drugage <- read_csv("./data/raw/drug_age/drugage.csv")
DA_drugs_incSyn <- unique(c(drugage$compound_name, 
                            drugage$synonyms))
DA_drugs_incSyn_CIDs <- convertIDs(DA_drugs_incSyn)
DA_drugs_noCID <- DA_drugs_incSyn_CIDs[[2]]
write_tsv(as.data.frame(DA_drugs_noCID), "./data/processed/drugAge_noCID.tsv")
DA_drugs_incSyn_CIDlist <- DA_drugs_incSyn_CIDs[[1]]
length(unique(DA_drugs_incSyn_CIDlist$drugName))
# [1] 369
length(unique(DA_drugs_incSyn_CIDlist$CID))
# [1] 526
nrow(unique(DA_drugs_incSyn_CIDlist))
# [1] 572
DA_drugs_incSyn_CHEMBLlist <- cid2chembl(DA_drugs_incSyn_CIDlist, 
                                         "./data/processed/drugAgeCHEMBL.tsv")
length(unique(DA_drugs_incSyn_CHEMBLlist$drugName))
# [1] 369
length(unique(DA_drugs_incSyn_CHEMBLlist$CID))
# [1] 526
length(setdiff(unique(DA_drugs_incSyn_CHEMBLlist$ChEMBLID), 
               NA))
# [1] 346
nrow(unique(DA_drugs_incSyn_CHEMBLlist))
# [1] 572

############### Adj Matrix

study_drug <- readRDS("./data/processed/study_drugMat.rds")
DA_drugs_incSyn_CHEMBLlist <- read_tsv("./data/processed/drugAgeCHEMBL.tsv")
drugage <- unique(setdiff(DA_drugs_incSyn_CHEMBLlist$ChEMBLID, 
                          NA))
drugs <- unique(c(colnames(study_drug), drugage))
drugage_mat <- rbind(rep(0, length(drugs)), 
                     rep(0, length(drugs)), as.numeric(drugs %in% drugage))
colnames(drugage_mat) <- drugs
rownames(drugage_mat) <- c("Human", "Model", "DrugAge")

genage_model1 <- read_tsv("./data/raw/genage_model/genage_models_orthologs_export.tsv") %>% 
    rename(model_symbol = `Model Organism Symbol`)
genage_model2 <- read_csv("./data/raw/genage_model/genage_models.csv") %>% 
    rename(model_symbol = symbol)
genage_model <- left_join(genage_model1, genage_model2, 
                          by = "model_symbol")
rm(genage_model1, genage_model2)
genage_humans <- read_csv("./data/raw/genage_human/genage_human.csv")
interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                        ChEMBLID) %>% unique()

drugs4genageHuman <- unique(setdiff((interactions %>% 
                                         filter(gene_name %in% genage_humans$symbol))$ChEMBLID, 
                                    NA))
drugs4genageModel <- unique(setdiff((interactions %>% 
                                         filter(gene_name %in% genage_model$Symbol))$ChEMBLID, 
                                    NA))

drugage_mat["Human", ] <- as.numeric(colnames(drugage_mat) %in% 
                                         drugs4genageHuman)
drugage_mat["Model", ] <- as.numeric(colnames(drugage_mat) %in% 
                                         drugs4genageModel)

saveRDS(drugage_mat, "./data/processed/lit_drugMat.rds")

################# Adj Mat for Genes
genage_model1 <- read_tsv("./data/raw/genage_model/genage_models_orthologs_export.tsv") %>% 
    rename(model_symbol = `Model Organism Symbol`)
genage_model2 <- read_csv("./data/raw/genage_model/genage_models.csv") %>% 
    rename(model_symbol = symbol)
genage_model <- left_join(genage_model1, genage_model2, 
                          by = "model_symbol")
rm(genage_model1, genage_model2)
genage_humans <- read_csv("./data/raw/genage_human/genage_human.csv")

study_gene <- readRDS("./data/processed/study_geneMat.rds")
genes <- unique(c(genage_humans$symbol, colnames(study_gene), 
                  genage_model$Symbol))

genage_genes <- rbind(as.numeric(genes %in% genage_humans$symbol), 
                      as.numeric(genes %in% genage_model$Symbol))

colnames(genage_genes) <- genes
rownames(genage_genes) <- c("Human", "Model")

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                        ChEMBLID) %>% unique()

DA_drugs_incSyn_CHEMBLlist <- read_tsv("./data/processed/drugAgeCHEMBL.tsv")
drugage <- unique(setdiff(DA_drugs_incSyn_CHEMBLlist$ChEMBLID, 
                          NA))
genes4drugage <- setdiff(unique(filter(interactions, 
                                       ChEMBLID %in% drugage)$gene_name), NA)

genage_genes <- rbind(genage_genes, as.numeric(colnames(genage_genes) %in% 
                                                   genes4drugage))
rownames(genage_genes)[3] <- "DrugAge"
genage_genes[, 1:10]
saveRDS(genage_genes, "./data/processed/lit_geneMat.rds")

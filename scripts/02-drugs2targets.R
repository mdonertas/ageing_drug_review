library(tidyverse)

mergeFiles <- function() {
    fx <- grep(".tsv|.rds|.csv", list.files("./data/processed/"), 
               invert = T, value = T)
    sapply(fx, function(fxx) {
        read_tsv(paste("./data/processed/", fxx, "/drugList_CHEMBL.tsv", 
                       sep = ""))
    })
}

allx <- mergeFiles()

noCHEMBL <- sapply(allx, function(x) {
    if ("drugName" %in% colnames(x)) {
        nms <- setdiff(unique(x$drugName[is.na(x$ChEMBLID)]), 
                       unique(x$drugName[!is.na(x$ChEMBLID)]))
        if (length(nms) == 0) {
            return("AllMapped")
        } else {
            nms
        }
    } else {
        return(NA)
    }
})

noCHEMBL <- noCHEMBL[noCHEMBL != "AllMapped" & !is.na(noCHEMBL)]
saveRDS(noCHEMBL, "./data/processed/noCHEMBL.rds")

cannotmap <- reshape2::melt(noCHEMBL) %>% rename(drugName = value, 
                                                 Study = L1) %>% unique()
cannotmap %>% group_by(Study) %>% summarise(n = length(unique(drugName)))
# A tibble: 3 x 2 Study n <chr> <int> 1 Barardo2017
# 5 2 Donertas2018 1 3 Liu2016 110
write_tsv(cannotmap, "./data/processed/noCHEMBL.tsv")

idList <- reshape2::melt(mergeFiles(), id.vars = c("ChEMBLID")) %>% 
    rename(other_id_type = variable, other_id = value, 
           study = L1) %>% select(ChEMBLID, study) %>% 
    unique() %>% na.omit()
write_tsv(idList, "./data/processed/mergedChemblIDs.tsv")
chemblIDs <- unique(idList$ChEMBLID)
length(chemblIDs)
# [1] 278

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% 
    select(gene_name, ChEMBLID, drug_name, interaction_claim_source) %>% 
    unique()

# do we have a match in dgiDB for the any of the
# drugs we couldn't map?
reshape2::melt(noCHEMBL) %>% rename(drug_name = value, 
                                    study = L1) %>% filter(study != "Liu2016") %>% 
    add_row(drug_name = "m40", study = "Barardo2017") %>% 
    mutate(drug_name = toupper(drug_name)) %>% left_join(interactions)
# no

combinedList <- left_join(idList, interactions) %>% 
    filter(!is.na(gene_name)) %>% unique()

write_tsv(combinedList, "./data/processed/combinedList.tsv")

combinedList %>% group_by(study) %>% 
    summarise(numDrugs = length(unique(ChEMBLID)),
              numGenes = length(unique(gene_name))) %>% knitr::kable()

# |study | numDrugs| numGenes|
# |:--------------|--------:|--------:| |Aliper2016
# | 7| 123| |Barardo2017 | 11| 113| |Calvert2016 |
# 7| 102| |Donertas2018 | 17| 261| |Fernandes2016 |
# 14| 92| |Fuentealba2018 | 10| 164| |Liu2016 | 54|
# 246| |Mofidifar2018 | 4| 5| |Snell2016 | 21| 150|
# |Snell2018 | 23| 140| |Yang2018 | 2| 8|
# |Ziehm2017 | 11| 103|

length(unique(combinedList$gene_name))
# 796

sum(chemblIDs %in% (combinedList %>% select(1, 3) %>% 
                        na.omit() %>% unique())[[1]])
# [1] 163

mean(chemblIDs %in% (combinedList %>% select(1, 3) %>% 
                         na.omit() %>% unique())[[1]])
# [1] 0.5863309

### EDA
drugOverlap <- sapply(unique(idList$study), function(x) {
    sapply(unique(idList$study), function(y) {
        x <- filter(idList, study == x)$ChEMBLID
        y <- filter(idList, study == y)$ChEMBLID
        length(unique(intersect(x, y)))/length(unique(x))
    })
})

diag(drugOverlap) <- NA
pheatmap::pheatmap(drugOverlap, display_numbers = T, 
                   breaks = seq(0, 1, length.out = 100), 
                   filename = "./results/drugOverlap.pdf", 
                   cellwidth = 20, cellheight = 20)

geneOverlap <- sapply(unique(combinedList$study), function(x) {
    sapply(unique(combinedList$study), function(y) {
        x <- unique(filter(combinedList, study == x)$gene_name)
        y <- unique(filter(combinedList, study == y)$gene_name)
        length(unique(intersect(x, y)))/length(unique(x))
    })
})

diag(geneOverlap) <- NA
pheatmap::pheatmap(geneOverlap, display_numbers = T, 
                   breaks = seq(0, 1, length.out = 100), 
                   filename = "./results/geneOverlap.pdf", 
                   cellwidth = 20, cellheight = 20)
### /EDA

### Construct Adj Matrices

study_drug <- combinedList %>% select(ChEMBLID, study) %>% 
    mutate(value = 1) %>% unique() %>% spread(ChEMBLID, 
                                              value, fill = 0)

gene_drug <- combinedList %>% select(ChEMBLID, gene_name) %>% 
    mutate(value = 1) %>% unique() %>% spread(ChEMBLID, 
                                              value, fill = 0)

study_gene <- combinedList %>% select(gene_name, study) %>% 
    mutate(value = 1) %>% unique() %>% spread(gene_name, 
                                              value, fill = 0)

sdmat <- as.matrix(select(study_drug, -study))
rownames(sdmat) <- study_drug$study

gdmat <- as.matrix(select(gene_drug, -gene_name))
rownames(gdmat) <- gene_drug$gene_name

sgmat <- as.matrix(select(study_gene, -study))
rownames(sgmat) <- study_gene$study

write_tsv(study_drug, "./data/processed/study_drugMat.tsv")
write_tsv(gene_drug, "./data/processed/gene_drugMat.tsv")
write_tsv(study_gene, "./data/processed/study_geneMat.tsv")

saveRDS(sdmat, "./data/processed/study_drugMat.rds")
saveRDS(gdmat, "./data/processed/gene_drugMat.rds")
saveRDS(sgmat, "./data/processed/study_geneMat.rds")

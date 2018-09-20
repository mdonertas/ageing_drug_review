# libraries
library(KEGGREST)
library(tidyverse)

## generate pathway data
kegg <- keggLink("pathway", "hsa")
kegg <- tibble(pathway = gsub("path:", "", kegg), 
               entrezgene = gsub("hsa:", "", names(kegg)))
martx <- biomaRt::useMart("ensembl", "hsapiens_gene_ensembl")
idmap <- biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene"), 
                        filters = "entrezgene", 
                        values = unique(kegg$entrezgene), mart = martx) %>% 
    mutate(entrezgene = as.character(entrezgene))

kegg <- left_join(kegg, idmap)
head(kegg)
rm(idmap)
kegg <- kegg %>% select(pathway, hgnc_symbol) %>% na.omit() %>% unique()
write_tsv(kegg, "./data/processed/kegg2gene.tsv")

study_gene <- readRDS("./data/processed/study_geneMat.rds")

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, ChEMBLID) %>% 
    unique()

study_pathway <- sapply(unique(kegg$pathway), function(pth) {
    allgns <- intersect(unique(filter(kegg, pathway == pth)$hgnc_symbol), 
                        interactions$gene_name)
    gns <- intersect(allgns, colnames(study_gene))
    if (length(gns) > 1) {
        mypath <- rowSums(study_gene[, gns])/length(allgns)
    } else if (length(gns) == 1) {
        mypath <- study_gene[, gns]/length(allgns)
    } else if (length(allgns) == 0) {
        mypath <- rep(NA, nrow(study_gene))
        names(mypath) <- rownames(study_gene)
    } else {
        mypath <- rep(0, nrow(study_gene))
        names(mypath) <- rownames(study_gene)
    }
    return(mypath)
})
names(which(colMeans(is.na(study_pathway)) == 1))
# [1] 'hsa00440' 'hsa00472' 'hsa00510' 'hsa00512' 'hsa00514' 'hsa00533'
# [7] 'hsa00534' 'hsa00601' 'hsa00780' 'hsa00785' 'hsa03020'
lit_genes <- readRDS("./data/processed/lit_geneMat.rds")

lit_pathway <- sapply(unique(kegg$pathway), function(pth) {
    allgns <- intersect(unique(filter(kegg, pathway == pth)$hgnc_symbol), 
                        interactions$gene_name)
    gns <- intersect(allgns, colnames(lit_genes))
    if (length(gns) > 1) {
        mypath <- rowSums(lit_genes[, gns])/length(allgns)
    } else if (length(gns) == 1) {
        mypath <- lit_genes[, gns]/length(allgns)
    } else if (length(allgns) == 0) {
        mypath <- rep(NA, nrow(lit_genes))
        names(mypath) <- rownames(lit_genes)
    } else {
        mypath <- rep(0, nrow(lit_genes))
        names(mypath) <- rownames(lit_genes)
    }
    return(mypath)
})
names(which(colMeans(is.na(lit_pathway)) == 1))
# [1] 'hsa00440' 'hsa00472' 'hsa00510' 'hsa00512' 'hsa00514' 'hsa00533'
# [7] 'hsa00534' 'hsa00601' 'hsa00780' 'hsa00785' 'hsa03020'
exc <- names(which(colMeans(is.na(lit_pathway)) == 1))
length(exc)
# 11
study_pathway <- study_pathway[, !colnames(study_pathway) %in% exc]
lit_pathway <- lit_pathway[, !colnames(lit_pathway) %in% exc]
dim(lit_pathway)
dim(study_pathway)
saveRDS(lit_pathway, "./data/processed/lit_pathwayMat.rds")
saveRDS(study_pathway, "./data/processed/study_pathwayMat.rds")

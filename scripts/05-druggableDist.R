library(scatterpie)
library(RFRlib)
library(ggthemes)
library(tidyverse)
library(ggrepel)
study_gene <- readRDS("./data/processed/study_geneMat.rds")
study_drug <- readRDS("./data/processed/study_drugMat.rds")
study_pathway <- readRDS("./data/processed/study_pathwayMat.rds")
martx <- biomaRt::useMart("ensembl", "hsapiens_gene_ensembl")

allgenes <- biomaRt::getBM(attributes = c("hgnc_symbol", "gene_biotype"), 
                           mart = martx)
unique(allgenes$gene_biotype)
allgenes <- filter(allgenes, gene_biotype == "protein_coding")$hgnc_symbol
length(allgenes)
# 19351
interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
    rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, ChEMBLID) %>% 
    unique()

length(unique(interactions$gene_name))
# 2994
geneswithinteractions <- unique(interactions$gene_name)

eachstsep_genes <- apply(study_gene, 1, function(x) {
    length(names(which(x == 1)))/length(geneswithinteractions)
})

numdisc <- colSums(study_gene)
ingenage <- names(numdisc) %in% read_csv("./data/raw/genage_human/genage_human.csv")$symbol
othergenes <- setdiff(geneswithinteractions, names(numdisc))
ingenage_others <- othergenes %in% read_csv("./data/raw/genage_human/genage_human.csv")$symbol

mydat <- data.frame(genes = c(names(numdisc), othergenes), 
                    numDisc = c(unname(numdisc), rep(0, length(othergenes))), 
                    inGenAge = c(ingenage, ingenage_others))

mydat <- mydat %>% unique() %>% group_by(numDisc, inGenAge) %>% 
    summarise(n = length(unique(genes))) %>% 
    mutate(inGenAge = c("Novel Discovery", "in GenAge")[1 + inGenAge]) %>% 
    ungroup()
mydat2 <- mydat %>% spread(inGenAge, n)
mydat2[is.na(mydat2)] <- 0
mydat2$n <- mydat2$`in GenAge` + mydat2$`Novel Discovery`

p1 <- ggplot()+
    geom_bar(data=mutate(mydat2,numDisc=numDisc*700),
             aes(x=numDisc,y=n),
             stat='identity',
             fill='gray20',width = 300)+
    geom_scatterpie(aes(x=numDisc,y=n+160),
                    data=mutate(mydat2,numDisc=numDisc*700),
                    cols = c('Novel Discovery','in GenAge'))+
    geom_label(data=mutate(mydat2,numDisc=numDisc*700),
               aes(x=numDisc+300,y=n+200,label=n))+
    scale_x_continuous(breaks=c(0:10)*700,labels = c(0:10))+
    theme_rfr()+
    scale_fill_manual(values=c('#FF9900','#3366CC'))+
    xlab('Number of Studies')+
    ylab('Number of Genes')+
    theme(legend.position = 'right')+
    guides(fill=guide_legend(''))+
    coord_fixed()

ggsave('./results/druggableGenomeDist.pdf',p1,width = 10,height = 5)


sort(numdisc,dec=T)[1:20]%>%
    knitr::kable()

# |       |  x|
# |:------|--:|
# |ABCB1  | 10|
# |BIRC5  |  9|
# |KRAS   |  9|
# |DDIT3  |  8|
# |ERBB2  |  8|
# |APC    |  7|
# |BAX    |  7|
# |BCL2   |  7|
# |BDNF   |  7|
# |BRAF   |  7|
# |EGFR   |  7|
# |FLT3   |  7|
# |MYC    |  7|
# |NTRK1  |  7|
# |PIK3CA |  7|
# |PIK3CG |  7|
# |RB1    |  7|
# |TP53   |  7|
# |AKT1   |  6|
# |CYP3A4 |  6|

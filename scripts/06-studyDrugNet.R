library(tidyverse)
druglist <- read_tsv("./data/processed/combinedList.tsv")
study_drug <- readRDS("./data/processed/study_drugMat.rds")
drugAge <- read_tsv("./data/processed/drugAgeCHEMBL.tsv")
aa <- filter(druglist, ChEMBLID %in% names(which(colSums(study_drug) > 1))) %>% 
    mutate(DrugAge = ChEMBLID %in% setdiff(unique(drugAge$ChEMBLID), NA)) %>% 
    select(drug_name, study, DrugAge) %>% unique() %>% arrange(drug_name)
library(igraph)
aa$drug_name <- tolower(aa$drug_name)
indrugAge <- unique(data.frame(aa$drug_name, aa$DrugAge))
indrugAge <- setNames(indrugAge[, 2], indrugAge[, 1])
myg <- graph_from_data_frame(aa)

library(GGally)
library(ggrepel)
V(myg)$type <- c("drug", "study")[1 + (V(myg)$name %in% unique(aa$study))]
stlab <- V(myg)$name
stlab[!stlab %in% unique(aa$study)] <- NA
drlab <- V(myg)$name
drlab[drlab %in% unique(aa$study)] <- NA
V(myg)$drugAge <- c("Novel", "in DrugAge")[1 + indrugAge[V(myg)$name]]
V(myg)$drugAge[is.na(V(myg)$drugAge)] <- "Study"
ggnet2(myg, size = 0, edge.color = "gray65", edge.size = 1, edge.alpha = 0.5, 
       layout.exp = 0) + 
    geom_point(size = 6, color = c("#3366CC", "#FF9900", NA)[factor(V(myg)$drugAge, 
                                                                    levels = c("Novel", "in DrugAge", "Study"))]) + 
    geom_label(label = stlab, fill = c(NA, "gray90")[factor(V(myg)$type)], 
               size = c(3, 5)[factor(V(myg)$type)]) + 
    geom_text_repel(label = drlab, color = "gray10", box.padding = 0.5, 
                    segment.size = 0)

ggsave("./results/st_drug_net.pdf", useDingbats = F, width = 10, height = 9)

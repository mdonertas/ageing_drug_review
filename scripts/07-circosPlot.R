## libraries
library(tidyverse)
library(KEGGREST)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

study_gene <- readRDS("./data/processed/study_geneMat.rds")
study_drug <- readRDS("./data/processed/study_drugMat.rds")
study_pathway <- readRDS("./data/processed/study_pathwayMat.rds")
genage_genes <- readRDS("./data/processed/lit_geneMat.rds")
drugage_mat <- readRDS("./data/processed/lit_drugMat.rds")
genage_pathway <- readRDS("./data/processed/lit_pathwayMat.rds")

subs <- colnames(study_gene) %in% (names(which(colSums(genage_genes[-3, 
                                                                    ]) == 0)))
study_gene[, subs][study_gene[, subs] == 1] <- 3
study_gene[study_gene == 1] <- 2
subs <- colnames(study_drug) %in% (names(which((drugage_mat[3, ]) == 0)))
study_drug[, subs][study_drug[, subs] == 1] <- 3
study_drug[study_drug == 1] <- 2

mat <- cbind(study_gene, study_drug, study_pathway)
rownames(mat) <- paste(rownames(mat), " ", sep = "")
mat <- mat[hclust(dist(mat))$order, ]

col_fun1 <- colorRamp2(c(seq(0, 1, length.out = 10), 2, 3), c("white", 
                                                              colorRampPalette(c("gray90", "gray25"))(9), "#FF9900", "#3366CC"))
lgd_pathway <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), col_fun = col_fun1, 
                      title = "KEGG Pathway Track")

col_fun2 <- colorRamp2(c(seq(0, 1, length.out = 10), 2), c("white", colorRampPalette(c("gray90", 
                                                                                       "gray25"))(9), "#DC3912"))

lgd_novel <- Legend(at = c("Novel Discovery"), type = "points", legend_gp = gpar(col = "#3366CC"), 
                    title_position = "topleft", title = "")

lgd_others <- Legend(at = c("Previously Discovered"), type = "points", 
                     legend_gp = gpar(col = "#FF9900"), title_position = "topleft", title = "")


lgd_list_vertical1 <- packLegend(lgd_pathway, lgd_novel, lgd_others)
lgd_list_vertical2 <- packLegend(lgd_pathway)
lgd_list_vertical3 <- packLegend(lgd_novel, lgd_others)

# > ncol(study_gene) [1] 796 > ncol(study_drug) [1] 163 >
# ncol(study_pathway) [1] 319
pdf("./results/circularHeatmap.pdf", width = 10.5, height = 10, useDingbats = F)

factors <- factor(c(rep("HUMAN GENES (796)", ncol(study_gene)), rep("DRUGS (163)", 
                                                                    ncol(study_drug)), rep("KEGG PATHWAYS (319)", ncol(study_pathway))), 
                  levels = c("DRUGS (163)", "HUMAN GENES (796)", "KEGG PATHWAYS (319)"))

mat_list <- list(`HUMAN GENES (796)` = mat[, factors == "HUMAN GENES (796)"], 
                 `DRUGS (163)` = mat[, factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[, 
                                                                                              factors == "KEGG PATHWAYS (319)"])

dend_list <- list(`HUMAN GENES (796)` = as.dendrogram(hclust(dist(t(mat_list[["HUMAN GENES (796)"]])))), 
                  `DRUGS (163)` = as.dendrogram(hclust(dist(t(mat_list[["DRUGS (163)"]])))), 
                  `KEGG PATHWAYS (319)` = as.dendrogram(hclust(dist(t(mat_list[["KEGG PATHWAYS (319)"]])))))



mat_list <- list(`HUMAN GENES (796)` = mat[c("Snell2016 ", "Snell2018 ", 
                                             "Mofidifar2018 "), factors == "HUMAN GENES (796)"], `DRUGS (163)` = mat[c("Snell2016 ", 
                                                                                                                       "Snell2018 ", "Mofidifar2018 "), factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[c("Snell2016 ", 
                                                                                                                                                                                                                 "Snell2018 ", "Mofidifar2018 "), factors == "KEGG PATHWAYS (319)"])
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = c(5, 5, 30), start.degree = 90, 
           clock.wise = T)
circos.initialize(factors, xlim = cbind(c(0, 0, 0), table(factors)))
circos.track(ylim = c(0, 3), track.height = 0.12, bg.border = "gray40", 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun1(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                             CELL_META$sector.index, niceFacing = TRUE, cex = 2, font = 2)
             })
circos.text(rep(0, nrow(mat_list[[1]])), 0:(nrow(mat_list[[1]]) - 1) + 
                1, rev(rownames(mat_list[[1]])), facing = "downward", niceFacing = T, 
            sector.index = "DRUGS (163)", cex = 1, adj = 1)

mat_list <- list(`HUMAN GENES (796)` = mat[c("Fernandes2016 ", "Fuentealba2018 "), 
                                           factors == "HUMAN GENES (796)"], `DRUGS (163)` = mat[c("Fernandes2016 ", 
                                                                                                  "Fuentealba2018 "), factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[c("Fernandes2016 ", 
                                                                                                                                                                               "Fuentealba2018 "), factors == "KEGG PATHWAYS (319)"])

circos.track(ylim = c(0, 2), track.height = 0.08, bg.border = "gray40", 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun1(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
             })
circos.text(rep(0, nrow(mat_list[[1]])), 0:(nrow(mat_list[[1]]) - 1) + 
                1, rev(rownames(mat_list[[1]])), facing = "downward", niceFacing = T, 
            sector.index = "DRUGS (163)", cex = 1, adj = 1)

mat_list <- list(`HUMAN GENES (796)` = mat[c("Liu2016 ", "Barardo2017 "), 
                                           factors == "HUMAN GENES (796)"], `DRUGS (163)` = mat[c("Liu2016 ", 
                                                                                                  "Barardo2017 "), factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[c("Liu2016 ", 
                                                                                                                                                                            "Barardo2017 "), factors == "KEGG PATHWAYS (319)"])

circos.track(ylim = c(0, 2), track.height = 0.08, bg.border = "gray40", 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun1(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
             })
circos.text(rep(0, nrow(mat_list[[1]])), 0:(nrow(mat_list[[1]]) - 1) + 
                1, rev(rownames(mat_list[[1]])), facing = "downward", niceFacing = T, 
            sector.index = "DRUGS (163)", cex = 1, adj = 1)


mat_list <- list(`HUMAN GENES (796)` = mat[c("Calvert2016 ", "Donertas2018 ", 
                                             "Yang2018 "), factors == "HUMAN GENES (796)"], `DRUGS (163)` = mat[c("Calvert2016 ", 
                                                                                                                  "Donertas2018 ", "Yang2018 "), factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[c("Calvert2016 ", 
                                                                                                                                                                                                          "Donertas2018 ", "Yang2018 "), factors == "KEGG PATHWAYS (319)"])

circos.track(ylim = c(0, 3), track.height = 0.12, bg.border = "gray40", 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun1(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
             })
circos.text(rep(0, nrow(mat_list[[1]])), 0:(nrow(mat_list[[1]]) - 1) + 
                1, rev(rownames(mat_list[[1]])), facing = "downward", niceFacing = T, 
            sector.index = "DRUGS (163)", cex = 1, adj = 1)

mat_list <- list(`HUMAN GENES (796)` = mat[c("Aliper2016 ", "Ziehm2017 "), 
                                           factors == "HUMAN GENES (796)"], `DRUGS (163)` = mat[c("Aliper2016 ", 
                                                                                                  "Ziehm2017 "), factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = mat[c("Aliper2016 ", 
                                                                                                                                                                          "Ziehm2017 "), factors == "KEGG PATHWAYS (319)"])

circos.track(ylim = c(0, 2), track.height = 0.08, bg.border = "gray40", 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun1(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
             })
circos.text(rep(0, nrow(mat_list[[1]])), 0:(nrow(mat_list[[1]]) - 1) + 
                1, rev(rownames(mat_list[[1]])), facing = "downward", niceFacing = T, 
            sector.index = "DRUGS (163)", cex = 1, adj = 1)


genage_genes <- readRDS("./data/processed/lit_geneMat.rds")
drugage_mat <- readRDS("./data/processed/lit_drugMat.rds")
genage_pathway <- readRDS("./data/processed/lit_pathwayMat.rds")
genage_genes[genage_genes == 1] <- 2
drugage_mat[drugage_mat == 1] <- 2
annotmat <- cbind(genage_genes, drugage_mat, genage_pathway)
annotmat <- annotmat[, colnames(mat)]
rownames(annotmat) <- paste(rownames(annotmat), " ", sep = "")
rownames(annotmat) <- paste("GenAge ", rownames(annotmat)[1:2], sep = "")

mat_list <- list(`HUMAN GENES (796)` = annotmat[, factors == "HUMAN GENES (796)"], 
                 `DRUGS (163)` = annotmat[, factors == "DRUGS (163)"], `KEGG PATHWAYS (319)` = annotmat[, 
                                                                                                        factors == "KEGG PATHWAYS (319)"])

circos.track(ylim = c(0, 3), bg.border = "gray40", track.height = 0.1, 
             panel.fun = function(x, y) {
                 sector.index <- CELL_META$sector.index
                 m <- mat_list[[sector.index]]
                 dend <- dend_list[[sector.index]]
                 
                 m2 <- m[, order.dendrogram(dend)]
                 col_mat <- col_fun2(m2)
                 nr <- nrow(m2)
                 nc <- ncol(m2)
                 for (i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 1:nc, rep(nr - i + 1, 
                                                                      nc), border = NA, col = col_mat[i, ])
                 }
             })
circos.text(rep(0, nrow(annotmat)), 0:(nrow(annotmat) - 1) + 1, rev(rownames(annotmat)), 
            facing = "downward", niceFacing = T, sector.index = "DRUGS (163)", 
            cex = 1, adj = 1)

max_height <- max(sapply(dend_list, function(x) attr(x, "height")))
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.25, 
             panel.fun = function(x, y) {
                 
                 sector.index <- get.cell.meta.data("sector.index")
                 dend <- dend_list[[sector.index]]
                 circos.dendrogram(dend, max_height = max_height)
             })

pushViewport(viewport(x = unit(2, "mm"), y = unit(6, "mm"), width = grobWidth(lgd_list_vertical1), 
                      height = grobHeight(lgd_list_vertical1), just = c("left", "bottom")))
grid.draw(lgd_list_vertical1)
upViewport()

circos.clear()
dev.off()

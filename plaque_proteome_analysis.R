setwd("c:/Users/kyang2/OneDrive - St. Jude Children's Research Hospital/Protein Turnover/plaque proteome")


library(readxl)
library(tidyverse)
library(ComplexHeatmap)
library(matrixStats)
library(Biobase)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(corrplot)
library(EnhancedVolcano)
library(circlize)
library(writexl)

###### fucntions for statistical analysis in this script ########
col_z_score <- function(x) {
        a <- t((x-colMeans(x)[col(x)]))/colSds(as.matrix(x))
        return(as.data.frame(t(a)))
}
row_z_score <- function(x) {
        a <- t((x-rowMeans(x)))/rowSds(as.matrix(x))[col(x)]
        return(as.data.frame(t(a)))
}
# limma package function for p.value and FDR calculation
limma.test <- function(data,range,design,cont.matrix,merge = "merge",label = "label") {
        merge = merge
        matrix <- data[,range] %>% as.matrix()
        row.names(matrix) <- data[,merge]
        eset.test <- ExpressionSet(assayData = matrix)
        fit <- lmFit(eset.test, design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        DE_table <- topTable(fit2, n = nrow(matrix),adjust="BH")
        colnames(DE_table)[c(1,4,5)] <- c("log2FC",
                                          "P Value",  
                                          "BH FDR")
        DE_table$`log2FC-z` <- col_z_score(DE_table[,1] %>% as.matrix())[[1]]
        colnames(DE_table) <- paste(label,colnames(DE_table),sep = " ")
        DE_table[,merge] <- rownames(DE_table)
        return(DE_table[,c(1,4,5,7,8)])
}


######## reading data and processing ##############
plaque_proteome <- read_xlsx("id_uni_prot_quan_includeMDK.xlsx")
names(plaque_proteome) <- as.vector(plaque_proteome[1,])
plaque_proteome <- plaque_proteome[-1,]
plaque_proteome[,24:54] <- lapply(plaque_proteome[,24:54],as.numeric)

names(plaque_proteome)[40:42] <- c("3M Log2FC plq v.s. nonplq",
                                   "3M p-value plq v.s. nonplq",
                                   "3M FDR plq v.s. nonplq")
names(plaque_proteome)[43:45] <- c("8M Log2FC plq v.s. nonplq",
                                   "8M p-value plq v.s. nonplq",
                                   "8M FDR plq v.s. nonplq")
names(plaque_proteome)[46:48] <- c("12M Log2FC plq v.s. nonplq",
                                   "12M p-value plq v.s. nonplq",
                                   "12M FDR plq v.s. nonplq")
names(plaque_proteome)[49:51] <- c("8M Log2FC plq v.s. wt",
                                   "8M p-value plq v.s. wt",
                                   "8M FDR plq v.s. wt")
names(plaque_proteome)[52:54] <- c("12M Log2FC plq v.s. wt",
                                   "12M p-value plq v.s. wt",
                                   "12M FDR plq v.s. wt")

names(plaque_proteome)

plaque_proteome_log2 <- plaque_proteome[,c(1:4,24:39)] %>% as.data.frame()
plaque_proteome_log2[,5:20] <- log2(plaque_proteome_log2[,5:20]) 

design <- matrix(data = (c(rep(c(1,0),4),rep(c(0,1),4))),nrow = 8, ncol = 2)
rownames(design) <- colnames(plaque_proteome_log2)[9:16]
colnames(design) <- c("plq","nonplq")
cont.matrix <- makeContrasts(plqvsnonplq=plq-nonplq, levels=design)
        
plaque_proteome_log2_limma_test <- limma.test(plaque_proteome_log2,
                                              range = c(9:16),
                                              design = design,
                                              cont.matrix = cont.matrix,
                                              merge = "Protein Accession #",
                                              label = "8M and 12M plq v.s. nonplq")

plaque_proteome_log2 <- merge(plaque_proteome_log2, plaque_proteome_log2_limma_test,
                              by = "Protein Accession #")

write_xlsx(plaque_proteome_log2[,c(1:4,21:24)],"plaque_proteome.xlsx")
plaque_proteome_log2_8m_12m <- plaque_proteome_log2[which(plaque_proteome_log2$`8M and 12M plq v.s. nonplq log2FC` >2 
                                                          &plaque_proteome_log2$`8M and 12M plq v.s. nonplq BH FDR` <0.05),]

plaque_proteome_log2_8m_12m<-plaque_proteome_log2_8m_12m[,c(1:4,
                                                            6,8,5,7,
                                                            10,12,9,11,
                                                            14,16,13,15,
                                                            17:24)]

plaque_8m_12m_top <- plaque_proteome_log2_8m_12m$GN %>% unique() %>% as.data.frame()
write_delim(plaque_8m_12m_top,"plaque_8m_12m_top.txt",col_names = F, na = "")

tiff("plaque_volcano.tiff", height = 1500, width = 1500, units = "px", res = 150)
EnhancedVolcano(plaque_proteome_log2[,c(4,21,23)],lab = plaque_proteome_log2$GN,
                x = "8M and 12M plq v.s. nonplq log2FC",
                xlab = "log2FC",
                y = "8M and 12M plq v.s. nonplq BH FDR",
                ylab = "-logFDR",
                pCutoff = 0.05,FCcutoff = 2,
                drawConnectors = T, subtitle = "8M and 12M plq v.s. nonplq")
dev.off()

col_fun <- colorRamp2(c(0,10),c("white","green"))
plaque_anno <- HeatmapAnnotation(Group = factor(rep(c("nonplq","nonplq","plq","plq"),3),
                                            ordered = T,
                                            levels = c("nonplq","plq")),
                            Month = factor(c(rep("3M",4),rep("8M",4),rep("12M",4)),
                                           ordered = T,
                                           levels = c("3M","8M","12M")),
                            col = list(Group = setNames(rainbow(2),c("nonplq","plq")),
                                       Month = setNames(col_fun(c(3,6,9)),c("3M","8M","12M"))),
                            border = T,
                            annotation_name_side = "right"
                            )

tiff("plaque_heatmap.tiff", width = 600, height = 1000,units = "px",res = 100)
Heatmap(as.matrix(plaque_proteome_log2_8m_12m[,5:16]%>% row_z_score()),
        show_column_names = FALSE,
        row_labels = plaque_proteome_log2_8m_12m$GN,
        show_row_dend = FALSE,
        row_names_side = "left",
        heatmap_legend_param = list(title = "z-score"),
        top_annotation = plaque_anno,
        cluster_columns = FALSE)
dev.off()


plaque_proteome_8m_12m <- plaque_proteome_log2[c(-grep("HUMAN",
                                                       plaque_proteome_log2$`Protein Accession #`),
                                                 -grep("CON_",
                                                       plaque_proteome_log2$`Protein Accession #`)),]

write_delim(as.data.frame(plaque_proteome_8m_12m$GN),"plaque_proteome_8m.txt",col_names = F)

plaque_proteome_8m_GO <- enrichGO(plaque_proteome_8m$GN,`org.Mm.eg.db`,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.01)
plaque_proteome_8m_GO_result <- plaque_proteome_8m_GO@result
plaque_proteome_8m_GO_pair <- pairwise_termsim(plaque_proteome_8m_GO)

plaque_proteome_8m_GO_fea <- emapplot(plaque_proteome_8m_GO_pair, cex_label_category=.8, cex_line=.5,showCategory = 15) + coord_cartesian() +
        scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                              guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
ggsave(plaque_proteome_8m_GO_fea, file="plaque_proteome_8m_GO_fea.tiff", width = 900,
       height = 700, units = "px",scale=3)



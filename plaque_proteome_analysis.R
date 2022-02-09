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
library(VennDiagram)

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
plaque_proteome <- read_xlsx("id_uni_prot_quan_includeMDK.xlsx") %>% as.data.frame()
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




# limma test
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

plaque_proteome_log2 <- plaque_proteome_log2 %>% group_by(GN) %>%
        slice_max(order_by = abs(`8M and 12M plq v.s. nonplq log2FC-z`),
                  n = 1, with_ties = F)

# pca analysis #
PCA1 <- prcomp(plaque_proteome_log2[,9:16]  %>% t(), scale. = TRUE)
PCA2 <- as.data.frame(PCA1$x) 
percentage <- round(PCA1$sdev / sum(PCA1$sdev) * 100, 2)
percentage <- paste(colnames(PCA2), "(", as.character(percentage), "%", ")", sep="")
PCA2$phenotype <- sapply(strsplit(rownames(PCA2),split = "_"),"[[",2)
PCA2$age <- sapply(strsplit(rownames(PCA2),split = "_"),"[[",3)
PCA2$age <- factor(PCA2$age, ordered = T, level = c("3m","8m","12m"))
PCA2$region <- c(rep(c("plq","nonplq"),4))
PCA2$region <- factor(PCA2$region, ordered = T, level = c("nonplq","plq"))


jpeg("plaque_PCA.jpeg",width = 1400, height = 1000, units = "px",res = 400)
ggplot(PCA2[1:8,],aes(x=PC1,y=PC2, color = region)) +
        geom_point(size = 4) +
        theme_classic() +
        xlab(percentage[1]) +
        ylab(percentage[2]) +
        scale_color_discrete()
dev.off()

##### check the important proteins ######
if(FALSE) {
Apoa2 <- plaque_proteome_log2 %>%
        filter(GN == "Apoa2") %>%
        pivot_longer(cols = 5:20,
                     names_to = "sample_name",
                     values_to = "log2TMT")
Apoa2 <- cbind(PCA2[,9:11],Apoa2[5:12,9:10])

jpeg("Apoa2.jpeg",width = 1200, height = 1000, units = "px",res = 200)
ggplot(data = Apoa2,aes(x = `age`, y = log2TMT, fill = region)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, position = "identity") 
dev.off()

Ralyl <- plaque_proteome_log2 %>%
        filter(GN == "Ralyl") %>%
        pivot_longer(cols = 5:20,
                     names_to = "sample_name",
                     values_to = "log2TMT")
Ralyl <- cbind(PCA2[,9:11],Ralyl[5:12,9:10])

jpeg("Ralyl.jpeg",width = 1200, height = 1000, units = "px",res = 200)
ggplot(data = Ralyl,aes(x = `age`, y = log2TMT, fill = region)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, position = "identity") 
dev.off()
}

## DE analysis ##
jpeg("plaque_heatmap.jpeg", height = 1500, width = 800, units = "px", res = 150)
Heatmap(plaque_proteome_log2[,c(10,12,14,16,
                                9,11,13,15)] %>% row_z_score() %>% as.matrix(),
        show_column_names = T,show_column_dend = F,
        show_row_names = F,show_row_dend = F)
dev.off()


plaque_proteome_log2_top <- plaque_proteome_log2 %>%
        filter(abs(`8M and 12M plq v.s. nonplq log2FC-z`) >2.5&
                       `8M and 12M plq v.s. nonplq BH FDR` <0.05) 



plaque_proteome_log2_top <- plaque_proteome_log2_top[order(
        abs(plaque_proteome_log2_top$`8M and 12M plq v.s. nonplq log2FC-z`),
        decreasing = T),]

jpeg("plaque_volcano.jpeg", height = 2000, width = 1500, units = "px", res = 200)
EnhancedVolcano(plaque_proteome_log2[,c("GN",
                                        "8M and 12M plq v.s. nonplq log2FC-z",
                                        "8M and 12M plq v.s. nonplq BH FDR")],
                lab = plaque_proteome_log2$GN,
                x = "8M and 12M plq v.s. nonplq log2FC-z",
                xlab = "log2FC-z (8M and 12M plq v.s. nonplq)",
                y = "8M and 12M plq v.s. nonplq BH FDR",
                ylab = "-logFDR", ylim = c(0,6),
                selectLab = c("Apoa2","Ralyl",plaque_proteome_log2_top$GN[1:15]),
                pCutoff = 0.05,FCcutoff = 2.5,maxoverlapsConnectors = 200,
                drawConnectors = T, subtitle = "8M and 12M plq v.s. nonplq")
dev.off()

plaque_proteome_log2_top2 <- plaque_proteome_log2 %>%
        filter(`8M and 12M plq v.s. nonplq log2FC-z` >2.5&
                       `8M and 12M plq v.s. nonplq BH FDR` <0.05)
names(plaque_proteome_log2_top2)[1] <- "Uniprot"

turnover_DE_slower <- read.delim("../turnover_analysis/turnover_DE_slower.txt",
                                 header = T,sep = "")

venn.diagram(
        x = list(turnover_DE_slower$GN,
                 plaque_proteome_log2_top2$GN),
        category.names = c("slower turnover proteins" , "plaque-enriched proteins" ),
        height = 1000, width = 1200, units = "px", resolution = 100,
        filename = 'venn_diagramm_nat_com.png',
        scaled = TRUE,
        ext.text = TRUE,
        ext.line.lwd = 2,
        ext.dist = -0.15,
        ext.length = 0.9,
        ext.pos = -4,
        inverted = TRUE,
        cex = 2.5,
        cat.cex = 2.5,
        rotation.degree = 45
)


turnover_palque_merge <- merge(turnover_DE_slower[,c("Uniprot","GN",
                                                     "genotype.anova.BH.FDR",
                                                     "Relative.half.life.change....")],
                               plaque_proteome_log2_top2[,c("Uniprot","GN",
                                                            "8M and 12M plq v.s. nonplq BH FDR",
                                                            "8M and 12M plq v.s. nonplq log2FC-z")],
                               by = c("Uniprot","GN"))

turnover_palque_merge[1,2] <- "A-beta"

library(cluster)
col_rhc = colorRamp2(c(0,100),
                     c( "white", "red"))
col_z = colorRamp2(c(0,max(turnover_palque_merge$`8M and 12M plq v.s. nonplq log2FC-z`)),
                       c("white", "red"))

jpeg("turnover_plaque_overlap.jpeg", height = 800, width = 500, units = "px", res = 200)
Heatmap(turnover_palque_merge[,4],
        col = col_rhc, row_labels = turnover_palque_merge$GN, 
        show_row_dend = F, row_names_side = "left",
        show_column_dend = F,show_column_names = F,
        cluster_rows = diana(turnover_palque_merge[,c(4,6)]),
        heatmap_legend_param = list(title = "")) +
        Heatmap(turnover_palque_merge[,6],
                col = col_z,show_column_names = F,
                heatmap_legend_param = list(title = ""))
dev.off()



# jumpn analysis
plaque_de_jumpn <- read.delim("z:/ResearchHome/ClusterHome/kyang2/JUMPn/plaque_DE/intermediate/PPI/foreground.txt_Publication_table_FDRfiltered.txt")
plaque_de_jumpn[,c(3:5,7:9)] <- plaque_de_jumpn[,c(3:5,7:9)] %>% 
        as.matrix() %>% 
        as.numeric()
plaque_de_jumpn$Database <- 
        sapply(strsplit(plaque_de_jumpn$GeneSetName, split = "_"), "[[",1)
plaque_de_jumpn <- plaque_de_jumpn %>% 
        mutate(BH = 0-log10(BH))

plaque_de_jumpn$Pathway <-
        str_replace(plaque_de_jumpn$GeneSetName,
                    pattern = "REACTOME_",
                    replacement = "") %>%
        str_replace(pattern = "HALLMARK_",
                    replacement = "") %>% 
        str_replace(pattern = "REACTOME_",
                    replacement = "") %>%
        str_replace(pattern = "GOBP_",
                    replacement = "") %>% 
        str_replace(pattern = "GOMF_",
                    replacement = "") %>%
        str_replace(pattern = "GOCC_",
                    replacement = "") %>%
        str_replace(pattern = "KEGG_",
                    replacement = "") %>%
        str_replace_all(pattern = "_",replacement = " ")

plaque_de_jumpn$Pathway <- plaque_de_jumpn$Pathway %>% tolower() %>% str_to_title()

if(FALSE) {
plaque_de_jumpn <- plaque_de_jumpn %>% 
        group_by(GeneSet) %>% 
        slice_max(order_by = BH, n = 1, with_ties = F)
}

c1m1 <- plaque_de_jumpn %>% 
        filter(ClusterName == "C1.M1") %>%
        mutate(ClusterName = "M1")
names(c1m1)[c(1:9)] <- paste(names(c1m1)[c(1:9)],"M1", sep = "_")
c1m2 <- plaque_de_jumpn %>% 
        filter(ClusterName == "C1.M2") %>%
        mutate(ClusterName = "M2")
names(c1m2)[c(1:9)] <- paste(names(c1m2)[c(1:9)],"M2",sep = "_")

c1m1m2 <- full_join(c1m1,c1m2, by = c("Pathway","Database"))
c1m1m2 <-c1m1m2[,c(grep("Pathway",names(c1m1m2)),
                   grep("Database",names(c1m1m2)),
                   grep("ClusterName",names(c1m1m2)),
                   grep("BH",names(c1m1m2)),
                   grep("GeneSet_",names(c1m1m2)))] %>%
        filter(BH_M1 > 1.3|BH_M2 >1.3)
c1m1m2$name_length <- nchar(c1m1m2$Pathway)



c1m1m2_select <- c1m1m2[unique(c(grep("APOE",c1m1m2$GeneSet_M1),
                                 grep("APP",c1m1m2$GeneSet_M2))),] 
c1m1m2_select <- c1m1m2_select[order(c1m1m2_select$BH_M1,na.last = T,decreasing = T),]

FDR_GO <- colorRamp2(c(0,2,4),c("white", "steelblue3","steelblue4"))
Heatmap(c1m1m2_select[,5:6] %>% as.matrix(),
        show_row_names = F,
        cluster_columns = F, cluster_rows = F,
        col = FDR_GO)


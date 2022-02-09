setwd("c:/Users/kyang2/OneDrive - St. Jude Children's Research Hospital/Protein Turnover/ProTransAnalysis")

library(tidyverse)
library(readxl)
library(Biobase)
library(limma)
library(ggstatsplot)
library(ComplexHeatmap)
library(viridisLite)
library(viridis)
library(circlize)
library(matrixStats)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(AnnotationDbi)
library(enrichplot)
library(writexl)
library(reshape2)
library(cluster)
##################### funtions to be used in this script ##################
# fucntion to calculate column-wise and row-wise z
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

### main analysis start here ###
########### read mouse dataset and data processing ###########
mouse_proteome <- read_excel("5xFAD_mouse_wholeProteome_PengLab_at_StJude_v1.0.0.xlsx",sheet = 'resorted')
colnames(mouse_proteome)[4:ncol(mouse_proteome)] = paste(colnames(mouse_proteome)[4:ncol(mouse_proteome)],"protein",sep = "_")
mouse_proteome[mouse_proteome == 0] <- NA
mouse_proteome <- mouse_proteome[complete.cases(mouse_proteome[4:19]),]
colnames(mouse_proteome)[1] <- "GN"

mouse_transcriptome <- read_excel("5xFAD_mouse_transcriptome_PengLab_at_StJude_v1.0.0.xlsx",sheet = 'resorted')
colnames(mouse_transcriptome)[4:ncol(mouse_transcriptome)] = paste(colnames(mouse_transcriptome)[4:ncol(mouse_transcriptome)],"mrna",sep = "_")
mouse_transcriptome[mouse_transcriptome == 0] <- NA
mouse_transcriptome <- mouse_transcriptome[complete.cases(mouse_transcriptome[4:18]),]
colnames(mouse_transcriptome)[1] <- "GN"
# log2 transformation
mouse_proteome_log2 <- rapply(mouse_proteome, f = log2, classes = c("numeric", "integer"), how = "replace")
mouse_proteome_log2 <- as.data.frame(mouse_proteome_log2)
mouse_transcriptome_log2 <- rapply(mouse_transcriptome, f = log2, classes = c("numeric", "integer"), how = "replace")
mouse_transcriptome_log2 <- as.data.frame(mouse_transcriptome_log2)
# data processing of proteome
design_proteome = matrix(c(rep(1,6),rep(0,6),rep(0,6),rep(1,6)),
                ncol = 2, byrow = F,
                nrow = ncol(mouse_proteome_log2[,grep("6mo|12mo",colnames(mouse_proteome_log2))]))
colnames(design_proteome) = c("FAD","WT")
rownames(design_proteome) = colnames(mouse_proteome_log2[grep("6mo|12mo",
                                                              colnames(mouse_proteome_log2))])
cont.matrix_proteome <- makeContrasts(FADvsWT=FAD-WT, levels=design_proteome)
mouse_proteome_limma_test <- limma.test(data = mouse_proteome_log2,
                range = grep("6mo|12mo",colnames(mouse_proteome_log2)),
                design = design_proteome,
                cont.matrix = cont.matrix_proteome,
                merge = "Protein Accession",
                label = "5XFAD v.s. WT protein")
mouse_proteome_log2 <- merge(mouse_proteome_log2,mouse_proteome_limma_test, by = "Protein Accession")
mouse_proteome_log2[which(mouse_proteome_log2$Annotation == "Amyloid beta peptide (K.LVFFAEDVGSNK.G)"),"GN"] <- "App"
if(FALSE) {
mouse_proteome_log2_unique <- mouse_proteome_log2 %>% group_by(GN) %>% 
  slice_max(order_by = abs(`5XFAD v.s. WT protein log2FC-z`), n = 1, with_ties = F)
}

# data processing of transcriptome
design_transcriptome = matrix(c(rep(1,6),rep(0,5),rep(0,6),rep(1,5)),
                         ncol = 2, byrow = F,
                         nrow = ncol(mouse_transcriptome_log2[,grep("6mo|12mo",
                                                                    colnames(mouse_transcriptome_log2))]))
colnames(design_transcriptome) = c("FAD","WT")
rownames(design_transcriptome) = colnames(mouse_transcriptome_log2[grep("6mo|12mo",
                                                                        colnames(mouse_transcriptome_log2))])
cont.matrix_transcriptome <- makeContrasts(FADvsWT=FAD-WT, levels=design_transcriptome)

mouse_transcriptome_limma_test <- limma.test(data = mouse_transcriptome_log2,
                                        range = grep("6mo|12mo",colnames(mouse_transcriptome_log2)),
                                        design = design_transcriptome,
                                        cont.matrix = cont.matrix_transcriptome,
                                        merge = "Gene ID",
                                        label = "5XFAD v.s. WT mRNA")
mouse_transcriptome_log2 <- merge(mouse_transcriptome_log2,mouse_transcriptome_limma_test, by = "Gene ID")
if(FALSE) {
mouse_transcriptome_log2_unique <- mouse_transcriptome_log2 %>% group_by(GN) %>% 
  slice_max(order_by = abs(`5XFAD v.s. WT mRNA log2FC-z`), n = 1, with_ties = F)
}

############# density plot of log2FC and log2FC-z ##########
ggplot() +
  geom_density(data = mouse_proteome_log2,
               aes(x= `5XFAD v.s. WT protein log2FC`, y = ..scaled.., fill = "red"),
               alpha = 0.5) +
  geom_density(data = mouse_transcriptome_log2,
               aes(x= `5XFAD v.s. WT mRNA log2FC`, y = ..scaled.., fill = "blue"),
               alpha = 0.5) +
  scale_fill_discrete(labels = c("protein","mRNA" )) +
  labs(fill = "readout", y = "Scaled density", x = "5XFAD v.s. WT log2FC") +
  theme_classic()

jpeg("pro_trans_colz_density.jpeg",width = 1000, height = 600, units = "px",res = 200)
ggplot() +
  geom_density(data = mouse_proteome_log2,
               aes(x= `5XFAD v.s. WT protein log2FC-z`, y = ..scaled.., fill = "red"),
               alpha = 0.5,size = 0.01) +
  geom_density(data = mouse_transcriptome_log2,
               aes(x= `5XFAD v.s. WT mRNA log2FC-z`, y = ..scaled.., fill = "blue"),
               alpha = 0.5,size = 0.01) +
  scale_fill_discrete(labels = c("protein","mRNA" )) +
  labs(fill = "readout", y = "Scaled density", x = "5XFAD v.s. WT log2FC-z") +
  theme_classic()
dev.off()

jpeg("pro_trans_colz_density_zoom.jpeg",width = 1000, height = 600, units = "px",res = 200)
ggplot() +
  geom_density(data = mouse_proteome_log2,
               aes(x= `5XFAD v.s. WT protein log2FC-z`, y = ..scaled.., fill = "red"),
               alpha = 0.5,size = 0.01) +
  geom_density(data = mouse_transcriptome_log2,
               aes(x= `5XFAD v.s. WT mRNA log2FC-z`, y = ..scaled.., fill = "blue"),
               alpha = 0.5,size = 0.01) +
  scale_fill_discrete(labels = c("protein","mRNA" )) +
  labs(fill = "readout", y = "Scaled density", x = "5XFAD v.s. WT log2FC-z") +
  theme_classic() +
  xlim(c(-2,2))
dev.off()

####### merge proteome and transcriptome data ##########
#pro_trans <- merge(mouse_proteome_log2_unique, mouse_transcriptome_log2_unique, by = "GN")
pro_trans <- full_join(mouse_proteome_log2, mouse_transcriptome_log2, by = "GN")
pro_trans <- pro_trans[!is.na(pro_trans$`Protein Accession`),]
pro_trans <- pro_trans[!is.na(pro_trans$`Gene ID`),]
names(pro_trans)
write_xlsx(pro_trans[,c(1:3,20:25,41:44)],"pro_trans.xlsx")


############ read human dataset and merge ########################
human_proteome <- read_excel("2021_ADProteomics_Review_Supp_v2.0.0.xlsx",sheet = 'resorted')
names(human_proteome)[1] <-"GN"
human_proteome$log2fc_z_protein <- col_z_score(human_proteome[,7] %>% as.matrix())[[1]]
human_proteome_distinct <- distinct(human_proteome,GN,.keep_all = TRUE)

human_transcriptome <- read_excel("MayoRNAseq_RNAseq_TCX_CONvsAD_DEG_Simple.xlsx")
human_transcriptome$beta_z_mrna <- col_z_score(human_transcriptome[,11] %>% as.matrix())[[1]]
names(human_transcriptome)[2] <- c("GN")
human_transcriptome_distinct <- distinct(human_transcriptome,GN,.keep_all = TRUE)

pro_trans_human <- merge(human_proteome_distinct, human_transcriptome_distinct, by = "GN")

write_xlsx(pro_trans_human[,c(1:3,7,6,8:9,10,20,19,22:24)],"pro_trans_human.xlsx")


######### mouse proteome transcriptome selection by z and FDR  ##########
# define the consistency by log2FC-z
for(i in 1:nrow(pro_trans)){
  if((pro_trans$`5XFAD v.s. WT mRNA log2FC-z`[i])*(pro_trans$`5XFAD v.s. WT protein log2FC-z`[i]) > 0) {
    if(abs(pro_trans$`5XFAD v.s. WT protein log2FC-z`[i] - pro_trans$`5XFAD v.s. WT mRNA log2FC-z`[i]) > 3
    ) {pro_trans$`consistency`[i] <- "inconsistent"
    }else{pro_trans$`consistency`[i] <- "consistent"}
  } else {
    pro_trans$`consistency`[i] <- "inconsistent"
  }}

col_z_protein = colorRamp2(c(min(pro_trans$`5XFAD v.s. WT protein log2FC-z`)
                             ,0,
                             max(pro_trans$`5XFAD v.s. WT protein log2FC-z`)),
                           c("blue", "white", "red"))
col_z_mrna = colorRamp2(c(min(pro_trans$`5XFAD v.s. WT mRNA log2FC-z`)
                          ,0,
                          max(pro_trans$`5XFAD v.s. WT mRNA log2FC-z`)),
                        c("blue", "white", "red"))

jpeg("z_heatmap_all.jpeg", height = 1000, width = 300, units = "px", res = 100)
Heatmap(pro_trans[order(pro_trans$`5XFAD v.s. WT protein log2FC-z`,
                               decreasing = T),
                  "5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F,cluster_rows = F,
        heatmap_legend_param = list(title = "Log2FC-z"),
        col = col_z_protein, border = T) +
  Heatmap(pro_trans[order(pro_trans$`5XFAD v.s. WT protein log2FC-z`,
                                 decreasing = T),
                    "5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = F,
          col = col_z_protein, border = T) +
  Heatmap(pro_trans[order(pro_trans$`5XFAD v.s. WT protein log2FC-z`,
                          decreasing = T),
                    "consistency"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          heatmap_legend_param = list(title = "consistency",
                                      border = T),
          col = c("white","black"),border = T)
dev.off()


pro_trans_top_up <- slice_max(filter(pro_trans,
                                            `5XFAD v.s. WT protein BH FDR` <0.05,
                                            `5XFAD v.s. WT protein log2FC-z` > 0),
                                     order_by = `5XFAD v.s. WT protein log2FC-z`,
                                     prop = 1)

pro_trans_top_down <- slice_min(filter(pro_trans,
                                              `5XFAD v.s. WT protein BH FDR` <0.05,
                                              `5XFAD v.s. WT protein log2FC-z` < 0),
                                       order_by = `5XFAD v.s. WT protein log2FC-z`,
                                       prop = 1)

if(FALSE) {
  consistentcy_count <- c(
    sum(pro_trans_top_up$consistency == "consistent")/nrow(pro_trans_top_up),
    1-sum(pro_trans_top_up$consistency == "consistent")/nrow(pro_trans_top_up),
    sum(pro_trans_top_down$consistency == "consistent")/nrow(pro_trans_top_down),
    1-sum(pro_trans_top_down$consistency == "consistent")/nrow(pro_trans_top_down)
  )
  print(consistentcy_count)
}

jpeg("z_heatmap_top_up.jpeg", height = 800, width = 200, units = "px", res = 100)
Heatmap(pro_trans_top_up[,"5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F,cluster_rows = F,
        show_heatmap_legend = F,
        col = col_z_protein, border = T,
        row_split =pro_trans_top_up[,"consistency"]) +
  Heatmap(pro_trans_top_up[,"5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = F,
          col = col_z_protein, border = T) 
dev.off()

jpeg("z_heatmap_top_down.jpeg", height = 300, width = 300, units = "px", res = 150)
Heatmap(pro_trans_top_down[,"5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_top_down$GN,
        row_names_side = "left",show_row_dend = F,
        show_heatmap_legend = F,
        col = col_z_protein, border = T) +
  Heatmap(pro_trans_top_down[,"5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = F,
          col = col_z_protein, border = T) +
  Heatmap(pro_trans_top_down[,"consistency"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = F,
          col = c("white","black"),border = T)
dev.off()

pro_trans_top_up_incon <- filter(pro_trans_top_up,
                                 consistency == "inconsistent") 
pro_trans_top_up_con <- filter(pro_trans_top_up,
                               consistency == "consistent") 
pro_trans_top_down_con <- filter(pro_trans_top_down,
                               consistency == "consistent")

pro_trans_top_up_incon$`z diff` <- pro_trans_top_up_incon$`5XFAD v.s. WT protein log2FC-z` - pro_trans_top_up_incon$`5XFAD v.s. WT mRNA log2FC-z`

pro_trans_top_up_incon_top <- pro_trans_top_up_incon%>%
  slice_max(order_by = `z diff`,prop = 0.5)

jpeg("z_heatmap_top_up_incon.jpeg", height = 1200, width = 250, units = "px", res = 150)
Heatmap(pro_trans_top_up_incon_top[,"5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_top_up_incon_top$GN,
        row_names_side = "left",show_row_dend = F,
        cluster_rows = pro_trans_top_up_incon_top[,c("5XFAD v.s. WT protein log2FC-z",
                                                     "5XFAD v.s. WT mRNA log2FC-z")] %>% diana(),
        show_heatmap_legend = F,
        col = col_z_protein, border = T) +
  Heatmap(pro_trans_top_up_incon_top[,"5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = F,
          col = col_z_protein, border = T) 
dev.off()

write.table(c(pro_trans_top_up_incon_top$`Protein Accession`,
              pro_trans_top_down$`Protein Accession`),
            "incon_proteins.txt",quote = F,row.names = F,
            col.names = T,)

enrich_GO_result <- function(a, group = "group") {
  a <- enrichGO(a,'org.Mm.eg.db',
                keyType = "SYMBOL",ont="ALL",
                pvalueCutoff=1,pAdjustMethod = "BH")
  a <- a@result
  names(a)[4:10] <- paste(group,names(a)[4:10], sep = "_")
  return(a)
}


pro_trans_top_up_incon_GO_result <- enrich_GO_result(pro_trans_top_up_incon$GN,
                                                     group = "up_inconsistent")
pro_trans_top_up_con_GO_result <- enrich_GO_result(pro_trans_top_up_con$GN,
                                                     group = "up_consistent")
pro_trans_top_down_con_GO_result <- enrich_GO_result(pro_trans_top_down_con$GN,
                                                     group = "down_consistent")

pro_trans_top_GO_result <- merge(pro_trans_top_up_incon_GO_result,
                                 pro_trans_top_up_con_GO_result,
                                 by = c("ID","Description","ONTOLOGY"),
                                 all = T)
pro_trans_top_GO_result <- merge(pro_trans_top_GO_result,
                                 pro_trans_top_down_con_GO_result,
                                 by = c("ID","Description","ONTOLOGY"),
                                 all = T)

pro_trans_top_GO_result <- pro_trans_top_GO_result[order(pro_trans_top_GO_result$up_inconsistent_p.adjust),c(1:3,
                                                      grep("p.adjust",names(pro_trans_top_GO_result)))]
pro_trans_top_GO_result[4:6] <- log10(pro_trans_top_GO_result[4:6])*(-1)
pro_trans_top_GO_result[is.na(pro_trans_top_GO_result)] <- 0

pro_trans_top_GO_result %>% 
  filter(ONTOLOGY == "MF") %>%
  slice_max(order_by = up_consistent_p.adjust,n = 10) 

GO_select <- c("lysosome",
               "transport vesicle",
               "early endosome",
               "synaptic vesicle",
               "lipoprotein particle",
               "amyloid precursor protein metabolic process",
               "regulation of amyloid-beta clearance",
               "tumor necrosis factor production",
               "response to lipoprotein particle",
               "amyloid-beta formation",
               "learning or memory",
               "regulation of endocytosis",
               "amyloid-beta binding",
               "heparin binding",
               "apolipoprotein binding",
               "MHC class I peptide loading complex",
               "extracellular membrane-bounded organelle",
               "endocytic vesicle",
               "response to interferon-gamma",
               "antigen processing and presentation",
               "cell killing",
               "TAP binding",
               "GTPase activity",
               "T cell receptor binding",
               "peptide binding")
               

FDR_GO <- colorRamp2(c(0,3,6),c("white", "steelblue3","steelblue4"))
 
pro_trans_top_GO_result_select <- filter(pro_trans_top_GO_result, 
                                         Description %in% GO_select)

jpeg("GO_heatmap.jpeg", width = 900, height = 900, units = "px", res = 150)
Heatmap(pro_trans_top_GO_result_select[,4:5] %>% 
          as.matrix(),
        row_labels = pro_trans_top_GO_result_select$Description,
        row_names_side = "right", show_column_names = F,
        cluster_rows = diana(pro_trans_top_GO_result_select[,4:5]),
        cluster_columns = F, show_row_dend = F,
        #row_split = pro_trans_top_GO_result_select$ONTOLOGY,
        col = FDR_GO, border = T,row_names_max_width = unit(12, "cm"),
        row_names_gp = gpar(fontsize = 16),
        heatmap_legend_param = list(title = "-LogFDR"))
dev.off()





# GO enrichment of top up inconsistent genes
if(FALSE) {
pro_trans_top_up_incon_GO <- enrichGO(pro_trans_top_up_incon$GN,'org.Mm.eg.db',
                                      keyType = "SYMBOL",ont="ALL",
                                      pvalueCutoff=0.05,pAdjustMethod = "BH")

pro_trans_top_up_incon_GO_result <- pro_trans_top_up_incon_GO@result %>% 
  filter(p.adjust < 0.001)


pro_trans_top_up_incon_GO<- pairwise_termsim(pro_trans_top_up_incon_GO)

pro_trans_top_up_incon_GO_amyloid <- grep("amyloid",pro_trans_top_up_incon_GO_result$Description)
pro_trans_top_up_incon_GO_amyloid <- c(pro_trans_top_up_incon_GO_amyloid,
                                       grep("memory",pro_trans_top_up_incon_GO_result$Description))
pro_trans_top_up_incon_GO_amyloid <- c(pro_trans_top_up_incon_GO_amyloid,
                                       grep("lipoprotein",pro_trans_top_up_incon_GO_result$Description))
pro_trans_top_up_incon_GO_amyloid <- unique(pro_trans_top_up_incon_GO_amyloid)
pro_trans_top_up_incon_GO_amyloid <- pro_trans_top_up_incon_GO_result[pro_trans_top_up_incon_GO_amyloid,
                                                                      "Description"]

jpeg("pro_trans_top_up_incon_GO.jpeg",width = 1000, height = 1400, units = "px",res = 150)
dotplot(pro_trans_top_up_incon_GO,
        showCategory = pro_trans_top_up_incon_GO_amyloid)+
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  theme(panel.grid.major.y = element_line(linetype='dotted',
                                          color='#808080'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  guides(size = guide_legend(override.aes=list(shape=1)))
dev.off()

jpeg("pro_trans_top_up_incon_GO_top20.jpeg",width = 1000, height = 1400, units = "px",res = 150)
dotplot(pro_trans_top_up_incon_GO,
        showCategory = 20)+
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  theme(panel.grid.major.y = element_line(linetype='dotted',
                                          color='#808080'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  guides(size = guide_legend(override.aes=list(shape=1)))
dev.off()
}

#### correlation of proteome and transcriptome #########

jpeg("pro_trans_colz_cor.jpeg",width = 800, height = 800, units = "px",res = 250)
ggscatterstats(
  data = pro_trans,
  x = "5XFAD v.s. WT mRNA log2FC-z",
  y = "5XFAD v.s. WT protein log2FC-z",
  xlab = "mRNA log2FC(5xFAD/WT)-z",
  ylab = "protein log2FC(5xFAD/WT)-z",
  marginal = F,
  results.subtitle = F,
  plotgrid.args = list(nrow = 3, ncol = 1),
  point.args = list(size = 2, alpha = 0.2, stroke = 0, colour = "blue"),
  smooth.line.args = list(size = 0.2, alpha = 0.5,
                          color = "pink", method = lm),
  bf.message = FALSE,
  ggplot.component = list(xlim(c(-15,25)), ylim(c(-15,25)))
)
dev.off()


jpeg("discrepancy.jpeg",width = 800, height = 800, units = "px",res = 300)
ggplot(pro_trans) +
  geom_point(aes(x = `5XFAD v.s. WT mRNA log2FC-z`,
                 y = `5XFAD v.s. WT protein log2FC-z`), 
             alpha = 1, size = 0.5) +
  xlim(c(-10,10))+
  ylim(c(-10,10)) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 1, alpha = 0.7)+
  geom_vline(xintercept = 0, size = 1, alpha = 0.7) +
  annotate("rect", xmin = -10, xmax = 0,
           ymin = 0, ymax = 10,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = 0, xmax = 10,
           ymin = -10, ymax = 0,
           alpha = .1,fill = "blue") +
  annotate("rect", xmin = 0, xmax = 10,
           ymin = 0, ymax = 10,
           alpha = .1,fill = "black") +
  annotate("rect", xmin = -10, xmax = 0,
           ymin = -10, ymax = 0,
           alpha = .1,fill = "black")
dev.off()

pdf("discrepancy.pdf",width = 4, height = 4)
ggplot(pro_trans) +
  geom_point(aes(x = `5XFAD v.s. WT mRNA log2FC-z`,
                 y = `5XFAD v.s. WT protein log2FC-z`), 
             alpha = 1, size = 0.5) +
  xlim(c(-10,10))+
  ylim(c(-10,10)) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 1, alpha = 0.7)+
  geom_vline(xintercept = 0, size = 1, alpha = 0.7) +
  annotate("rect", xmin = -10, xmax = 0,
           ymin = 0, ymax = 10,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = 0, xmax = 10,
           ymin = -10, ymax = 0,
           alpha = .1,fill = "blue") +
  annotate("rect", xmin = 0, xmax = 10,
           ymin = 0, ymax = 10,
           alpha = .1,fill = "black") +
  annotate("rect", xmin = -10, xmax = 0,
           ymin = -10, ymax = 0,
           alpha = .1,fill = "black")
dev.off()

jpeg("pro_trans_human_colz_cor.jpeg",width = 800, height = 800, units = "px",res = 150)
ggscatterstats(
  data = pro_trans_human,
  x = "beta_z_mrna",
  y = "log2fc_z_protein",
  label.var = GN,
  label.expression = (GN == "APP" |
                        GN == "APOE"|
                        GN == "SLIT2"|
                        GN == "SFRP1"|
                        GN == "SMOC1"|
                        GN == "HTRA1"|
                        GN == "MDK"|
                        GN == "NTN1"|
                        GN == "CTHRC1"|
                        GN == "NTN3"|
                        GN == "SLT3"|
                        GN == "LSP1"|
                        GN == "C4B"|
                        GN == "CLU"|
                        GN == "OLFML3"|
                        GN == "ICAM1"),
  xlab = "mRNA beta-z",
  ylab = "protein log2FC(AD/ctl)-z",
  results.subtitle = F,
  marginal = F,
  point.args = list(size = 2, alpha = 0.2, stroke = 0,
                    colour = "blue", method = lm),
  smooth.line.args = list(size = 0.2, alpha = 0.5,
                          color = "pink", method = lm),
  bf.message = FALSE,
  ggplot.component = list(xlim(c(-10,25)), ylim(c(-10,25)))
)
dev.off()


if(FALSE) {
######### mouse proteome transcriptome selection by FDR and z  ##########
pro_trans_scaled <- pro_trans[,c(1:3,22:23,43:44)]
names(pro_trans_scaled)
pro_trans_scaled$`5XFAD v.s. WT protein -log(FDR)` <- pro_trans_scaled$`5XFAD v.s. WT protein BH FDR` %>% log10() %>% abs()
pro_trans_scaled$`5XFAD v.s. WT mRNA -log(FDR)` <- pro_trans_scaled$`5XFAD v.s. WT mRNA BH FDR` %>% log10() %>% abs()

jpeg("FDR_distribution_all.jpeg",height = 800, width = 1000, units = "px", res = 250)
ggplot(melt(pro_trans_scaled[,8:9],variable.name = "readout", value.name = "-log(FDR)")) +
  geom_density(aes(x = `-log(FDR)`, y = ..scaled.., color = `readout`)) +
  scale_color_manual(labels = c("protein", "mRNA"), values = c("blue", "red"))
dev.off()

jpeg("FDR_distribution_top.jpeg",height = 800, width = 1000, units = "px", res = 250)
ggplot(melt(pro_trans_scaled[,8:9],variable.name = "readout", value.name = "-log(FDR)")) +
  geom_density(aes(x = `-log(FDR)`, y = ..scaled.., color = `readout`)) +
  scale_color_manual(labels = c("protein", "mRNA"), values = c("blue", "red")) +
  xlim(c(0,2))
dev.off()

for(i in 1:nrow(pro_trans_scaled)){
  if(pro_trans_scaled$`5XFAD v.s. WT protein log2FC-z`[i] > 0) {
    pro_trans_scaled$`5XFAD v.s. WT protein adjust -log(FDR)`[i] <- 
      pro_trans_scaled$`5XFAD v.s. WT protein -log(FDR)`[i]*(1)
    } else {
      pro_trans_scaled$`5XFAD v.s. WT protein adjust -log(FDR)`[i] <-
        pro_trans_scaled$`5XFAD v.s. WT protein -log(FDR)`[i]*(-1)
    }}

for(i in 1:nrow(pro_trans_scaled)){
  if(pro_trans_scaled$`5XFAD v.s. WT mRNA log2FC-z`[i] > 0) {
    pro_trans_scaled$`5XFAD v.s. WT mRNA adjust -log(FDR)`[i] <- 
      pro_trans_scaled$`5XFAD v.s. WT mRNA -log(FDR)`[i]*(1)
  } else {
    pro_trans_scaled$`5XFAD v.s. WT mRNA adjust -log(FDR)`[i] <-
      pro_trans_scaled$`5XFAD v.s. WT mRNA -log(FDR)`[i]*(-1)
  }}

for(i in 1:nrow(pro_trans_scaled)){
  if((pro_trans_scaled$`5XFAD v.s. WT mRNA log2FC-z`[i])*(pro_trans_scaled$`5XFAD v.s. WT protein log2FC-z`[i]) > 0) {
    if(abs(pro_trans_scaled$`5XFAD v.s. WT protein log2FC-z`[i] - pro_trans_scaled$`5XFAD v.s. WT mRNA log2FC-z`[i]) > 3
    ) {pro_trans_scaled$`consistency`[i] <- "inconsistent"
    }else{pro_trans_scaled$`consistency`[i] <- "consistent"}
  } else {
    pro_trans_scaled$`consistency`[i] <- "inconsistent"
  }}

col_FDR_protein = colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_FDR_mrna = colorRamp2(c(-3, 0, 6), c("blue", "white", "red"))

jpeg("FDR_heatmap_all.jpeg", height = 1000, width = 200, units = "px", res = 100)
Heatmap(pro_trans_scaled[order(pro_trans_scaled$`5XFAD v.s. WT protein adjust -log(FDR)`,
                               decreasing = T),
                         10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F,cluster_rows = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T) +
  Heatmap(pro_trans_scaled[order(pro_trans_scaled$`5XFAD v.s. WT protein adjust -log(FDR)`,
                                 decreasing = T),
                           11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna, border = T)
dev.off()

pro_trans_scaled_top_up <- slice_max(filter(pro_trans_scaled,
                                            `5XFAD v.s. WT protein BH FDR` <0.05,
                                            `5XFAD v.s. WT protein log2FC-z` > 0),
                                     order_by = `5XFAD v.s. WT protein adjust -log(FDR)`,
                                     prop = 1)

pro_trans_scaled_top_down <- slice_min(filter(pro_trans_scaled,
                                              `5XFAD v.s. WT protein BH FDR` <0.05,
                                              `5XFAD v.s. WT protein log2FC-z` < 0),
                                       order_by = `5XFAD v.s. WT protein adjust -log(FDR)`,
                                       prop = 1)

if(FALSE) {
consistentcy_count <- c(
sum(pro_trans_scaled_top_up$consistency == "consistent")/nrow(pro_trans_scaled_top_up),
1-sum(pro_trans_scaled_top_up$consistency == "consistent")/nrow(pro_trans_scaled_top_up),
sum(pro_trans_scaled_top_down$consistency == "consistent")/nrow(pro_trans_scaled_top_down),
1-sum(pro_trans_scaled_top_down$consistency == "consistent")/nrow(pro_trans_scaled_top_down)
)
print(consistentcy_count)
}

jpeg("FDR_heatmap_top_up.jpeg", height = 600, width = 150, units = "px", res = 80)
Heatmap(pro_trans_scaled_top_up[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F, show_row_dend = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T,show_heatmap_legend = F, 
        row_split =pro_trans_scaled_top_up[,12])+
  Heatmap(pro_trans_scaled_top_up[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna,border = T,show_heatmap_legend = F) 
dev.off()

jpeg("FDR_heatmap_top_down.jpeg", height = 200, width = 200, units = "px", res = 100)
Heatmap(pro_trans_scaled_top_down[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_scaled_top_down$GN,
        row_names_side = "left",show_row_dend = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T,show_heatmap_legend = F,
        row_split =pro_trans_scaled_top_down[,12])+
  Heatmap(pro_trans_scaled_top_down[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna,border = T,show_heatmap_legend = F) 
dev.off()

pro_trans_scaled_top_up_incon <- filter(pro_trans_scaled_top_up,
                                        consistency == "inconsistent") %>%
  slice_max(order_by = `5XFAD v.s. WT protein adjust -log(FDR)`,prop = 0.5)

write_delim(pro_trans_scaled_top_up_incon$GN %>% as.data.frame(),
            "top_up_inconsistent_genes.txt",col_names = F)

pro_trans_scaled_top_up_con <- filter(pro_trans_scaled_top_up,
                                        consistency == "consistent") %>%
  slice_max(order_by = `5XFAD v.s. WT protein adjust -log(FDR)`,prop = 0.1)

jpeg("FDR_heatmap_top_up_inconsistent.jpeg", height = 900, width = 230, units = "px", res = 120)
Heatmap(pro_trans_scaled_top_up_incon[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_scaled_top_up_incon$GN,
        row_names_side = "left",show_row_dend = F,
        cluster_rows =diana(pro_trans_scaled_top_up_incon[,10:11]),
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T,show_heatmap_legend = F)+
  Heatmap(pro_trans_scaled_top_up_incon[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna,border = T,show_heatmap_legend = F) 
dev.off()

jpeg("FDR_heatmap_top_up_consistent.jpeg", height = 250, width = 200, units = "px", res = 110)
Heatmap(pro_trans_scaled_top_up_con[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_scaled_top_up_con$GN,
        row_names_side = "left",show_row_dend = F,
        cluster_rows =diana(pro_trans_scaled_top_up_con[,10:11]),
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T,show_heatmap_legend = F)+
  Heatmap(pro_trans_scaled_top_up_con[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna,border = T,show_heatmap_legend = F) 
dev.off()


######### human proteome transcriptome selection by FDR and z  ##########
pro_trans_human_scaled <- pro_trans_human[,c(1:3,8:9,23:24)]
pro_trans_human_scaled$`protein -log(FDR)` <- pro_trans_human_scaled$`BH FDR` %>% log10() %>% abs()
pro_trans_human_scaled$`mRNA -log(FDR)` <- pro_trans_human_scaled$Dx.qValue %>% log10() %>% abs()

jpeg("FDR_distribution_human_all.jpeg",height = 800, width = 1000, units = "px", res = 250)
ggplot(melt(pro_trans_human_scaled[,8:9],variable.name = "readout", value.name = "-log(FDR)")) +
  geom_density(aes(x = `-log(FDR)`, y = ..scaled.., color = `readout`)) +
  scale_color_manual(labels = c("protein", "mRNA"), values = c("blue", "red"))
dev.off()

jpeg("FDR_distribution_human_top.jpeg",height = 800, width = 1000, units = "px", res = 250)
ggplot(melt(pro_trans_human_scaled[,8:9],variable.name = "readout", value.name = "-log(FDR)")) +
  geom_density(aes(x = `-log(FDR)`, y = ..scaled.., color = `readout`)) +
  scale_color_manual(labels = c("protein", "mRNA"), values = c("blue", "red")) +
  xlim(c(0,7))
dev.off()

for(i in 1:nrow(pro_trans_human_scaled)){
  if(pro_trans_human_scaled$log2fc_z_protein[i] > 0) {
    pro_trans_human_scaled$`protein adjust -log(FDR)`[i] <- 
      pro_trans_human_scaled$`protein -log(FDR)`[i]*(1)
  } else {
    pro_trans_human_scaled$`protein adjust -log(FDR)`[i] <- 
      pro_trans_human_scaled$`protein -log(FDR)`[i]*(-1)
  }}

for(i in 1:nrow(pro_trans_human_scaled)){
  if(pro_trans_human_scaled$beta_z_mrna[i] > 0) {
    pro_trans_human_scaled$`mRNA adjust -log(FDR)`[i] <- 
      pro_trans_human_scaled$`mRNA -log(FDR)`[i]*(1)
  } else {
    pro_trans_human_scaled$`mRNA adjust -log(FDR)`[i] <- 
      pro_trans_human_scaled$`mRNA -log(FDR)`[i]*(-1)
  }}

for(i in 1:nrow(pro_trans_human_scaled)){
  if((pro_trans_human_scaled$beta_z_mrna[i])*(pro_trans_human_scaled$log2fc_z_protein[i]) > 0) {
    if(abs(pro_trans_human_scaled$log2fc_z_protein[i] - pro_trans_human_scaled$beta_z_mrna[i]) > 3
    ) {pro_trans_human_scaled$`consistency`[i] <- "inconsistent"
    }else{pro_trans_human_scaled$`consistency`[i] <- "consistent"}
  } else {
    pro_trans_human_scaled$`consistency`[i] <- "inconsistent"
  }}

col_FDR_human_protein = colorRamp2(c(-22, 0, 31), c("blue", "white", "red"))
col_FDR_human_mrna = colorRamp2(c(-9, 0, 9), c("blue", "white", "red"))

jpeg("FDR_heatmap_human_all.jpeg", height = 1000, width = 200, units = "px", res = 100)
Heatmap(pro_trans_human_scaled[order(pro_trans_human_scaled$`protein adjust -log(FDR)`,
                               decreasing = T),10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F,cluster_rows = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_human_protein, border = T) +
  Heatmap(pro_trans_human_scaled[order(pro_trans_human_scaled$`protein adjust -log(FDR)`,
                                 decreasing = T),11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna, border = T)
dev.off()

pro_trans_human_scaled_top_up <- slice_max(filter(pro_trans_human_scaled,
                                            `BH FDR` <0.05,
                                            log2fc_z_protein > 0),
                                     order_by = `protein adjust -log(FDR)`,
                                     prop = 1)

pro_trans_human_scaled_top_down <- slice_min(filter(pro_trans_human_scaled,
                                              `BH FDR` <0.05,
                                              log2fc_z_protein < 0),
                                       order_by = `protein adjust -log(FDR)`,
                                       prop = 1)

if(FALSE) {
  consistentcy_count <- c(
    sum(pro_trans_human_scaled_top_up$consistency == "consistent")/nrow(pro_trans_human_scaled_top_up),
    1-sum(pro_trans_human_scaled_top_up$consistency == "consistent")/nrow(pro_trans_human_scaled_top_up),
    sum(pro_trans_human_scaled_top_down$consistency == "consistent")/nrow(pro_trans_human_scaled_top_down),
    1-sum(pro_trans_human_scaled_top_down$consistency == "consistent")/nrow(pro_trans_human_scaled_top_down)
  )
  print(consistentcy_count)
}

jpeg("FDR_heatmap_human_top_up.jpeg", height = 600, width = 150, units = "px", res = 80)
Heatmap(pro_trans_human_scaled_top_up[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F, show_row_dend = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_human_protein, border = T,show_heatmap_legend = F, 
        row_split =pro_trans_human_scaled_top_up[,12])+
  Heatmap(pro_trans_human_scaled_top_up[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna,border = T,show_heatmap_legend = F) 
dev.off()

jpeg("FDR_heatmap_human_top_down.jpeg", height = 600, width = 150, units = "px", res = 80)
Heatmap(pro_trans_human_scaled_top_down[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        show_row_names = F,show_row_dend = F,
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_human_protein, border = T,show_heatmap_legend = F,
        row_split =pro_trans_human_scaled_top_down[,12])+
  Heatmap(pro_trans_human_scaled_top_down[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna,border = T,show_heatmap_legend = F) 
dev.off()

pro_trans_human_scaled_top_up_incon <- filter(pro_trans_human_scaled_top_up,
                                        consistency == "inconsistent") %>%
  slice_max(order_by = `protein adjust -log(FDR)`,prop = 0.25)
write_delim(pro_trans_human_scaled_top_up_incon$GN %>% as.data.frame(),
            "human_top_up_inconsistent_genes.txt",col_names = F)

pro_trans_human_scaled_top_down_incon <- filter(pro_trans_human_scaled_top_down,
                                              consistency == "inconsistent") %>%
  slice_min(order_by = `protein adjust -log(FDR)`,prop = 0.1)
write_delim(pro_trans_human_scaled_top_down_incon$GN %>% as.data.frame(),
            "human_top_down_inconsistent_genes.txt",col_names = F)


jpeg("FDR_heatmap_human_top_up_inconsistent.jpeg", height = 1100, width = 210, units = "px", res = 100)
Heatmap(pro_trans_human_scaled_top_up_incon[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_human_scaled_top_up_incon$GN,
        row_names_side = "left",show_row_dend = F,
        cluster_rows =diana(pro_trans_human_scaled_top_up_incon[,10:11]),
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_human_protein, border = T,show_heatmap_legend = F)+
  Heatmap(pro_trans_human_scaled_top_up_incon[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna,border = T,show_heatmap_legend = F) 
dev.off()

jpeg("FDR_heatmap_human_top_down_inconsistent.jpeg", height = 1100, width = 210, units = "px", res = 100)
Heatmap(pro_trans_human_scaled_top_down_incon[,10] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = pro_trans_human_scaled_top_down_incon$GN,
        row_names_side = "left",show_row_dend = F,
        cluster_rows =diana(pro_trans_human_scaled_top_down_incon[,10:11]),
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_human_protein, border = T,show_heatmap_legend = F)+
  Heatmap(pro_trans_human_scaled_top_down_incon[,11] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna,border = T,show_heatmap_legend = F) 
dev.off()

####### overlap of mouse and human inconsistent genes ############
human_maRt = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_maRt = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
huamn_mouse_GNs <- getLDS(attributes = "external_gene_name", 
                          filters = "external_gene_name", 
                          values = pro_trans_human_scaled_top_up$GN, 
                          mart = human_maRt, 
                          attributesL = "external_gene_name", 
                          martL = mouse_maRt) 
names(huamn_mouse_GNs) <- c("GN","Mouse GN")
pro_trans_human_scaled_top_up <- merge(pro_trans_human_scaled_top_up,
                                       huamn_mouse_GNs, by = "GN")
pro_trans_scaled_top_up$`Mouse GN` <- pro_trans_scaled_top_up$GN

mouse_human_top_up <- merge(pro_trans_scaled_top_up,
                            pro_trans_human_scaled_top_up,
                            by = "Mouse GN")
mouse_human_top_up <- filter(mouse_human_top_up,
                             consistency.x == "inconsistent",
                             consistency.y == "inconsistent",)

jpeg("FDR_heatmap_mouse_human_top_up_inconsistent.jpeg",
     height = 200, width = 200, units = "px", res = 80)
Heatmap(mouse_human_top_up[,11] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = mouse_human_top_up$`Mouse GN`,
        row_names_side = "left",show_row_dend = F,
        cluster_rows =diana(mouse_human_top_up[,c(11:12,23:24)]),
        heatmap_legend_param = list(title = "-Log(FDR)"),
        col = col_FDR_protein, border = T,show_heatmap_legend = F)+
  Heatmap(mouse_human_top_up[,12] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_mrna,border = T,show_heatmap_legend = F) +
  Heatmap(mouse_human_top_up[,23] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_protein, border = T,show_heatmap_legend = F)+
  Heatmap(mouse_human_top_up[,24] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F, show_row_dend = F,
          heatmap_legend_param = list(title = "-Log(FDR)"),
          col = col_FDR_human_mrna,border = T,show_heatmap_legend = F) 
dev.off()

}
################ protein-based selection ###########################
protein_up <- pro_trans[which(pro_trans$`5XFAD v.s. WT protein log2FC-z` >2
                                  &pro_trans$`5XFAD v.s. WT mRNA log2FC-z` <1
                                  &pro_trans$`5XFAD v.s. WT protein BH FDR` < 0.05
),]

write_delim(as.data.frame(protein_up[,1]),"protein_up.txt",col_names = F)

protein_up_GO <- enrichGO(protein_up[,1],'org.Mm.eg.db', keyType = "SYMBOL",ont="BP",
                          pvalueCutoff=0.05)

protein_up_GO_result <- protein_up_GO@result


protein_up_GO<- pairwise_termsim(protein_up_GO)

jpeg("protein_up_GO_fea.jpeg",width = 1000, height = 800, units = "px",res = 150)
emapplot(protein_up_GO,repel = T, #node_label ="none",
         cex_label_category=.8, cex_line=0.5,
         showCategory = 15) + 
  coord_cartesian() +
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                        guide = guide_colorbar(reverse = TRUE, order=1),
                        trans='log10')
dev.off()

write.csv(protein_up_GO_result[1:15,],"protein_up_GO_topGNs.csv")
protein_up_GO_topGNs <- strsplit(protein_up_GO_result[1:15,8]
                                 ,split = "/") %>% unlist() %>% unique()
print(protein_up_GO_topGNs)
print(length(protein_up_GO_topGNs))

jpeg("protein_up_cor.jpeg",width = 800, height = 800, units = "px",res = 200)
ggscatterstats(
  # arguments relevant for ggstatsplot::ggscatterstats
  data = pro_trans,
  x = "5XFAD v.s. WT mRNA log2FC-z",
  y = "5XFAD v.s. WT protein log2FC-z",
  label.var = GN,
  label.expression = c(pro_trans$GN %in% protein_up_GO_topGNs),
  xlab = "mRNA log2FC(5xFAD/WT)-z",
  ylab = "protein log2FC(5xFAD/WT)-z",
  results.subtitle = F,
  marginal = F,
  plotgrid.args = list(nrow = 3, ncol = 1),
  point.args = list(size = 2, alpha = 0.2, stroke = 0, colour = "blue"),
  smooth.line.args = list(size = 0.2, alpha = 0.5,
                          color = "pink",method = lm),
  bf.message = FALSE,
  ggplot.component = list(ggplot2::annotate("rect", xmin = -10, xmax = 1,
                                            ymin = 2, ymax = 24,
                                            alpha = .1,fill = "red"),
  xlim(c(-10,15)), ylim(c(-10,25)))
)
dev.off()

jpeg("protein_down_cor.jpeg",width = 800, height = 800, units = "px",res = 150)
ggscatterstats(
  data = pro_trans,
  x = "5XFAD v.s. WT mRNA log2FC-z",
  y = "5XFAD v.s. WT protein log2FC-z",
  label.var = GN,
  label.expression = (`5XFAD v.s. WT protein log2FC-z` < -1
                        &`5XFAD v.s. WT mRNA log2FC-z` > -0.5
                        &`5XFAD v.s. WT protein BH FDR` < 0.05
  ),
  xlab = "mRNA log2FC(5xFAD/WT)-z",
  ylab = "protein log2FC(5xFAD/WT)-z",
  results.subtitle = F,
  marginal = F,
  plotgrid.args = list(nrow = 3, ncol = 1),
  point.args = list(size = 2, alpha = 0.2, stroke = 0, colour = "blue"),
  smooth.line.args = list(size = 0.2, alpha = 0.5,
                          color = "pink", method = lm),
  bf.message = FALSE,
  ggplot.component = list(ggplot2::annotate("rect", xmin = -0.5, xmax = 20, ymin = -10, ymax = -1,
                                            alpha = .1,fill = "red"),
                          xlim(c(-10,25)), ylim(c(-10,25)))
)
dev.off()

# heatmap of selected genes in row_z
protein_up_z <- cbind(row_z_score(protein_up[,4:19]),row_z_score(protein_up[,26:40]))

row.names(protein_up_z) <- protein_up[,1]
protein_up_z <- protein_up_z[(row.names(protein_up_z) %in% protein_up_GO_topGNs),
                             ] 
protein_up_z <- protein_up_z[,c(9:16,1:8,25:31,17:24)] # resort the data for heatmap

col_fun <- colorRamp2(c(-12,0,12), c("red", "white", "blue"))
pro_tran_anno <- HeatmapAnnotation(readout = c(rep("protein",16),rep("mRNA",15)) %>% factor(levels = c("protein","mRNA")),
                                   genotype = c(rep("WT",8),rep("5XFAD",8),rep("WT",7),rep("5XFAD",8)),
                                   month = c(rep("3-month",2),rep("6-month",4),rep("12-month",2),
                                             rep("3-month",2),rep("6-month",4),rep("12-month",2),
                                             rep("3-month",2),rep("6-month",3),rep("12-month",2),
                                             rep("3-month",2),rep("6-month",3),rep("12-month",3)
                                             ) %>% factor(levels = c("3-month","6-month","12-month")),                                            
                                   col = list(readout = setNames(col_fun(c(-12,12)),c("protein","mRNA")) ,
                                              genotype = setNames(col_fun(c(-6,6)),c("WT","5XFAD")),
                                              month = setNames(col_fun(c(3,6,12)),c("3-month","6-month","12-month"))),
                                   border = T,
                                   simple_anno_size = unit(0.3, "cm"),
                                   annotation_name_side = "right")
jpeg("protein_up_heatmap.jpeg",width = 1600, height = 1000, units = "px",res = 200)                             
Heatmap(protein_up_z,cluster_columns = F, column_names_gp = gpar(fontsize = 10),
        show_column_names = F,show_row_names = F,show_row_dend = F, show_heatmap_legend = F)
dev.off()

# heatmap of selected genes in log2fc-z
protein_up_top <- protein_up[protein_up$GN %in% protein_up_GO_topGNs,]
pro_tran_anno2 <- HeatmapAnnotation(readout = c("protein","mRNA"),
                                    col = list(readout = setNames(col_fun(c(-12,12)),c("protein","mRNA"))),
                                    border = T,
                                    simple_anno_size = unit(1, "cm"),
                                    annotation_name_side = "left")
col_fun2 = colorRamp2(c(-2, 0, 24), c("green", "white", "red"))

jpeg("protein_up_heatmap2.jpeg",width = 400, height = 600, units = "px",res = 100)
Heatmap(protein_up_top[,c(23,44)] %>%as.matrix(),
        cluster_columns = F, show_column_names = F,
        row_names_side = "left",row_labels = protein_up_top$GN,show_row_dend = F,
        heatmap_legend_param = list(title = "Log2FC z-score"),
        top_annotation = pro_tran_anno2,col = col_fun2
        )
dev.off()


# combiend heatmap     
jpeg("protein_up_heatmap3.jpeg",width = 1600, height = 800, units = "px",res = 200)
Heatmap(protein_up_top[,c(23,44)] %>%as.matrix(),
          cluster_columns = F, show_column_names = F,
          row_names_side = "left",row_labels = protein_up_top$GN,show_row_dend = F,
          heatmap_legend_param = list(title = "Log2FC z-score"),
          top_annotation = pro_tran_anno2,col = col_fun2
  )+
    Heatmap(protein_up_z,cluster_columns = F, column_names_gp = gpar(fontsize = 10),
        show_column_names = F,row_names_side = "left",show_row_dend = F,
        heatmap_legend_param = list(title = "TMT-intensity/FKPB z-score"),
        top_annotation = pro_tran_anno) 
dev.off()

protein_up_overlap <- read.delim("../turnover_analysis/overlap.txt")
protein_up_overlap <- protein_up_overlap[,1]
protein_up_overlap2 <- protein_up[protein_up$GN %in% protein_up_overlap,]
protein_up_z2 <- cbind(row_z_score(protein_up[,4:19]),row_z_score(protein_up[,26:40]))

row.names(protein_up_z2) <- protein_up[,1]
protein_up_z2 <- protein_up_z2[(row.names(protein_up_z2) %in% protein_up_overlap),
] 
protein_up_z2 <- protein_up_z2[,c(9:16,1:8,25:31,17:24)] # resort the data for heatmap

jpeg("protein_up_overlap_heatmap.jpeg",width = 1600, height = 800, units = "px",res = 200)
Heatmap(protein_up_overlap2[,c(23,44)] %>%as.matrix(),
        cluster_columns = F, show_column_names = F,
        row_names_side = "left",row_labels = protein_up_overlap2$GN,show_row_dend = F,
        heatmap_legend_param = list(title = "Log2FC z-score"),
        top_annotation = pro_tran_anno2,col = col_fun2
)+
  Heatmap(protein_up_z2,cluster_columns = F, column_names_gp = gpar(fontsize = 10),
          show_column_names = F,row_names_side = "left",show_row_dend = F,
          heatmap_legend_param = list(title = "TMT-intensity/FKPB z-score"),
          top_annotation = pro_tran_anno) 
dev.off()

# human dataset
protein_up_human <- pro_trans_human[which(pro_trans_human$log2fc_z_protein >2
                              &pro_trans_human$beta_z_mrna <1
                              &pro_trans_human$`BH FDR` < 0.05
                              #&pro_trans_human$p.value_6m_12m_mrna < 0.05
),]

protein_up_human_GO <- enrichGO(protein_up_human[,1],'org.Hs.eg.db', keyType = "SYMBOL",ont="BP",
                          pvalueCutoff=0.05)
protein_up_human_GO_result <- protein_up_human_GO@result
protein_up_human_GO_topGNs <- strsplit(protein_up_human_GO_result[1:30,8]
                                 ,split = "/") %>% unlist() %>% unique()
print(protein_up_human_GO_topGNs)
print(length(protein_up_human_GO_topGNs))

protein_up_human_GO <- pairwise_termsim(protein_up_human_GO)

jpeg("protein_up_human_GO_fea.jpeg",width = 1000, height = 800, units = "px",res = 125)
emapplot(protein_up_human_GO, cex_label_category=.8, cex_line=.5,showCategory = 30) + coord_cartesian() +
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
dev.off()

jpeg("protein_up_human_cor.jpeg",width = 800, height = 800, units = "px",res = 150)
ggscatterstats(
  # arguments relevant for ggstatsplot::ggscatterstats
  data = pro_trans_human,
  x = "beta_z_mrna",
  y = "log2fc_z_protein",
  label.var = GN,
  label.expression = (pro_trans_human$GN %in% protein_up_human_GO_topGNs),
  xlab = "mRNA beta-z",
  ylab = "protein log2FC(AD/ctl)-z",
  results.subtitle = F,
  marginal = F,
  point.args = list(size = 2, alpha = 0.2, stroke = 0,
                    colour = "blue", method = lm),
  smooth.line.args = list(size = 0.2, alpha = 0.5,
                          color = "pink",method = lm),
  bf.message = FALSE,
  ggplot.component = list(ggplot2::annotate("rect", xmin = -10, xmax = 1, ymin = 2, ymax = 17,
                                            alpha = .1,fill = "red"),
                          xlim(c(-10,25)), ylim(c(-10,25)))
)
dev.off()

################ mRNA-based selection ###########################
jpeg("mrna_up_cor.jpeg",width = 800, height = 800, units = "px",res = 150)
ggscatterstats(
  # arguments relevant for ggstatsplot::ggscatterstats
  data = pro_trans,
  x = "5XFAD v.s. WT mRNA log2FC-z",
  y = "5XFAD v.s. WT protein log2FC-z",
  label.var = GN,
  label.expression = (`5XFAD v.s. WT protein log2FC-z` < 0.5
                      &`5XFAD v.s. WT mRNA log2FC-z` > 2
                      &`5XFAD v.s. WT mRNA BH FDR` < 0.05
  ),
  xlab = "mRNA log2FC(5xFAD/WT)-z",
  ylab = "protein log2FC(5xFAD/WT)-z",
  results.subtitle = F,
  marginal = F,
  plotgrid.args = list(nrow = 3, ncol = 1),
  point.args = list(size = 2, alpha = 0.2, stroke = 0, colour = "blue"),
  smooth.line.args = list(size = 0.2, alpha = 0.5, color = "pink"),
  bf.message = FALSE,
  ggplot.component = list(ggplot2::annotate("rect", xmin = 2, xmax = 25, ymin = -10, ymax = 0.5,
                                            alpha = .1,fill = "red"),
                          xlim(c(-10,25)), ylim(c(-10,25)))
)
dev.off()

jpeg("mrna_down_cor.jpeg",width = 800, height = 800, units = "px",res = 150)
ggscatterstats(
  # arguments relevant for ggstatsplot::ggscatterstats
  data = pro_trans,
  x = "5XFAD v.s. WT mRNA log2FC-z",
  y = "5XFAD v.s. WT protein log2FC-z",
  label.var = GN,
  label.expression = (`5XFAD v.s. WT protein log2FC-z` > -0.5
                      &`5XFAD v.s. WT mRNA log2FC-z` < -1
                      &`5XFAD v.s. WT mRNA BH FDR` < 0.05
  ),
  xlab = "mRNA log2FC(5xFAD/WT)-z",
  ylab = "protein log2FC(5xFAD/WT)-z",
  results.subtitle = F,
  marginal = F,
  plotgrid.args = list(nrow = 3, ncol = 1),
  point.args = list(size = 2, alpha = 0.2, stroke = 0, colour = "blue"),
  smooth.line.args = list(size = 0.2, alpha = 0.5, color = "pink"),
  bf.message = FALSE,
  ggplot.component = list(ggplot2::annotate("rect", xmin = -10, xmax = -1, ymin = -0.5, ymax = 25,
                                            alpha = .1,fill = "red"),
                          xlim(c(-10,25)), ylim(c(-10,25)))
)
dev.off()

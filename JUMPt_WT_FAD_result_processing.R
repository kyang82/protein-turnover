setwd("c:/Users/kyang2/OneDrive - St. Jude Children's Research Hospital/Protein Turnover/turnover_analysis/")

library(readxl)
library(tidyverse)
library(ComplexHeatmap)
library(ggstatsplot)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(writexl)
library(colorspace)
library(reshape2)
library(circlize)
library(limma)
library(Biobase)
########## dataset processing and overview of turnover profile ##############
wt_fad_turnover <- read.delim("uni_protein_turnover_final.txt",
                              sep = "\t",header = T,skip = 1) %>% as.data.frame()
names(wt_fad_turnover)

wt_fad_turnover <- wt_fad_turnover[,c(21,1,20,3:17,23:37)]

names(wt_fad_turnover)[4:33] <- c(paste(rep("wt",10),c("day-0", "day-4", "day-8", "day-16", "day-32"),1:10, sep = "_"),
                                  paste(rep("fad",5),c("day-0", "day-4", "day-8", "day-16", "day-32"),1:5, sep = "_"),
                                  paste(rep("wt",5),c("day-0", "day-4", "day-8", "day-16", "day-32"),11:15, sep = "_"),
                                  paste(rep("fad",10),c("day-0", "day-4", "day-8", "day-16", "day-32"),6:15, sep = "_"))
names(wt_fad_turnover)

### case study of apoe ###
if(FALSE) {
wt_fad_turnover_melt <- pivot_longer(wt_fad_turnover[1116,],
                                     cols = 4:33,
                                     names_to = "sample",
                                     values_to = "L%",
                                     values_drop_na = T) %>% as.data.frame()
wt_fad_turnover_melt <- wt_fad_turnover_melt[which(wt_fad_turnover_melt$`L%` != 1),]

wt_fad_turnover_melt$genotype <- sapply(str_split(wt_fad_turnover_melt$sample,"_"),"[[",1)
wt_fad_turnover_melt$time <- factor(sapply(str_split(wt_fad_turnover_melt$sample,"_"),"[[",2),
                                    ordered = T, level = c("day-0", "day-4", "day-8", "day-16", "day-32"))
ggplot(wt_fad_turnover_melt) +
  geom_boxplot(aes(x = time, y = `L%`, fill = genotype)) 

res.aov2 <- aov(`L%` ~ genotype*time, data = wt_fad_turnover_melt)
summary(res.aov2)[[1]][1:2,5]
wt_fad_turnover$genotype_anova[1] <- summary(res.aov2)[[1]][1:2,5][1]
wt_fad_turnover$time_anova[1] <- summary(res.aov2)[[1]][1:2,5][2]
}

### function for two-way ANOVA  ###
two_way_anova <- function(a){
for ( i in 1:nrow(a)) {
melt <- pivot_longer(a[i,],
                     cols = 4:33,
                     names_to = "sample",
                     values_to = "L%",
                     values_drop_na = T) %>% as.data.frame()
melt <- melt[which(melt$`L%` != 1),]
melt$genotype <- sapply(str_split(melt$sample,"_"),"[[",1)
melt$time <- factor(sapply(str_split(melt$sample,"_"),"[[",2),
                    ordered = T,
                    level = c("day-0", "day-4", "day-8", "day-16", "day-32"))
res.aov2 <- aov(`L%` ~ genotype*time, data = melt)
a$`genotype anova p value`[i] <- summary(res.aov2)[[1]][1:2,5][1]
a$`time anova p value`[i] <- summary(res.aov2)[[1]][1:2,5][2]
}
  a$`genotype anova BH FDR` <- p.adjust(p = a$`genotype anova p value`,
                                        method = "BH")
  return(a)}

wt_fad_turnover <- two_way_anova(wt_fad_turnover)
sum(wt_fad_turnover$`genotype anova p value` < 0.05)
sum(wt_fad_turnover$`genotype anova BH FDR` < 0.05)

ggplot(wt_fad_turnover) +
  geom_density(aes(x = `genotype anova p value`, y = ..scaled..,
                   fill = "blue"), alpha = 0.5) +
  geom_density(aes(x = `genotype anova BH FDR`, y = ..scaled..,
                   fill = "red"), alpha = 0.5) +
  scale_fill_discrete(labels = c("p value","FDR" )) +
  labs(fill = "Two-way ANOVA", y = "Scaled density", x = "5XFAD v.s. WT log2FC") +
  geom_vline(xintercept = 0.05) +
  xlim(c(0,1))



names(wt_fad_turnover)[2] <- "Uniprot"

########## dataset processing and overview of JUMPt result ##############
### reading the file and sort the dataset ###
wt_fad_halflives <- read_xlsx("results_SILAC_TMT_JUMPt_11_27_w_cutoff_wo_lys.xlsx",sheet = "results") %>% as.data.frame()
A_beta <- wt_fad_halflives[18060,]
wt_fad_halflives <- wt_fad_halflives[-18060,]
wt_fad_halflives$Uniprot <- sapply(strsplit(wt_fad_halflives$Var1,'\\$'),"[[",1)
wt_fad_halflives$mice_type <- sapply(strsplit(wt_fad_halflives$Var1,'\\$'),"[[",2)
wt_fad_halflives$mice_type <- factor(wt_fad_halflives$mice_type, levels = c("WT","FAD"),ordered = TRUE)
names(wt_fad_halflives)[c(8,12)] <- c("Half-life (days)","mice type")
names(wt_fad_halflives)[2] <- "GN"
### boxplot to overview both half-lifes of WT and FAD ###
tiff("wt_fad_overview_boxplot.png",height = 800, width = 600,units = "px",res = 100)
ggplot(data = wt_fad_halflives) +
  geom_boxplot(aes(x= `mice type`,y = `Half-life (days)`,fill = `mice type`),show.legend = F) +
  xlab("Mice Type") +
  ylab("Protein Halflife (days)") 
dev.off()

tiff("wt_fad_overview_density.tiff",width = 1000, height = 600, units = "px",res = 200)
ggplot(wt_fad_halflives) +
  geom_density(aes(x= `Half-life (days)`,
                   y = ..scaled..,
                   fill = `mice type`),
               alpha = 0.5,size = 0.01) +
  labs(fill = "mice type", y = "Scaled density", x = "Half-life (days)") +
  theme_classic() +
  xlim(c(0,20))
dev.off()


### manage data to seperate and then merge WT amd FAD datasets ###
wt_halflives <- wt_fad_halflives[which(wt_fad_halflives$`mice type` == "WT"),]
wt_halflives <- wt_halflives[,c("Uniprot","GN","Half-life (days)",
                                "time_point_1","time_point_2","time_point_3","time_point_4","time_point_5")]
names(wt_halflives)[4:8] <- c("day-0", "day-4", "day-8", "day-16", "day-32")
fad_halflives <- wt_fad_halflives[which(wt_fad_halflives$`mice type` == "FAD"),]
fad_halflives <- fad_halflives[,c("Uniprot","GN","Half-life (days)",
                                "time_point_1","time_point_2","time_point_3","time_point_4","time_point_5")]
names(fad_halflives)[4:8] <- c("day-0", "day-4", "day-8", "day-16", "day-32")

wt_fad_merge <- merge(wt_halflives,fad_halflives,by = c("Uniprot","GN"))
names(wt_fad_merge)[3:8] <- paste("WT",c("","day-0", "day-4", "day-8", "day-16", "day-32"))
names(wt_fad_merge)[9:14] <- paste("FAD",c("","day-0", "day-4", "day-8", "day-16", "day-32"))
names(wt_fad_merge)[c(3,9)] <- c("WT","FAD")

wt_fad_merge_simple <- merge(wt_fad_merge[,c(1:3,9)],wt_fad_turnover[,c(2,34,35)],by = c("Uniprot"))

### calculate the log2 fold change of half-lives between FAD and WT mice ###
wt_fad_merge$`log2fc of halflife (FAD vs WT)` <- log2(wt_fad_merge$`FAD`/wt_fad_merge$`WT`)
log2fc <- wt_fad_merge$`log2fc of halflife (FAD vs WT)`
wt_fad_merge$`log2fc z` <- (log2fc-mean(log2fc))/sd(log2fc)
remove(wt_halflives,fad_halflives,log2fc)
wt_fad_merge$mean.halflife <- (wt_fad_merge$WT + wt_fad_merge$FAD )/2
wt_fad_merge$UNIPROT1 <- sapply(strsplit(wt_fad_merge$Uniprot, "\\|"),"[[",2)
wt_fad_merge$UNIPROT2 <- as.character(mapIds(org.Mm.eg.db,keys = wt_fad_merge$GN,
                                keytype = "SYMBOL",column = "UNIPROT"))
wt_fad_merge$GeneID1 <- mapIds(org.Mm.eg.db,keys = wt_fad_merge$UNIPROT1,
                           keytype = "UNIPROT",column = "ENTREZID")
wt_fad_merge$GeneID2 <- mapIds(org.Mm.eg.db,keys = wt_fad_merge$UNIPROT2,
                               keytype = "UNIPROT",column = "ENTREZID")

###### master table of four datasets #########
pro_trans <- read_xlsx("../ProTransAnalysis/pro_trans.xlsx")
pro_trans <- pro_trans[,c(1:3,7,5:6,13,11:12)]

plaque_proteome <- read_xlsx("../plaque proteome/plaque_proteome.xlsx")
names(plaque_proteome)[1] <- c("Protein Accession")
plaque_proteome <- plaque_proteome[,c(1,8,6:7)]

turnover <- wt_fad_merge[,c(1:3,9)]
names(turnover)[c(1,3,4)] <- c("Protein Accession", "WT Half-life (days)","5XFAD Half-life (days)")

turnover_meta_table <- merge(pro_trans,plaque_proteome, by = c("Protein Accession"), all =T)
turnover_meta_table <- merge(turnover_meta_table,turnover,by = c("Protein Accession","GN"), all =T)
turnover_meta_table <- turnover_meta_table[-grep("CON_",turnover_meta_table$`Protein Accession`),]
names(turnover_meta_table)[1:3] <- paste("Mouse",names(turnover_meta_table)[1:3], sep = " ")

pro_trans_human <- read_xlsx(("../ProTransAnalysis/pro_trans_human.xlsx"))

human_maRt = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_maRt = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
huamn_mouse_GNs <- getLDS(attributes = "external_gene_name", 
                             filters = "external_gene_name", 
                             values = pro_trans_human$GN, 
                             mart = human_maRt, 
                             attributesL = "external_gene_name", 
                             martL = mouse_maRt) 
names(huamn_mouse_GNs) <- c("GN","Mouse GN")
pro_trans_human <- merge(pro_trans_human,huamn_mouse_GNs, by = "GN")
names(pro_trans_human)[1:3] <- paste("Human",names(pro_trans_human)[1:3], sep = " ")
names(pro_trans_human)[4:7] <- paste("AD v.s. ctl protein", 
                                     c("log2FC", "P value", "BH FDR","log2FC-z"))
names(pro_trans_human)[10:13] <- paste("AD v.s. ctl mRNA", 
                                     c("log2FC", "P value", "BH FDR","log2FC-z"))
pro_trans_human <- pro_trans_human[,c(1:3,7,5:6,13,11:12,14)]

turnover_meta_table <- merge(pro_trans_human,turnover_meta_table, by = "Mouse GN")

write_xlsx(turnover_meta_table,"turnover_meta_table.xlsx")

########### overview of L% #############

# selected L% heatmap of both mice types
days = colorRamp2(c(4, 32), c("white", "blue"))
wt_fad_anno <- HeatmapAnnotation(Phenotype = factor(c("WT","WT","WT","WT","5XFAD","5XFAD","5XFAD","5XFAD"),
                                                    ordered = TRUE,levels = c("WT","5XFAD")),
                                 "Labeling time (days)" = c(4,8,16,32,4,8,16,32),
                                 col = list(Phenotype = setNames(rainbow(2),c("WT","5XFAD")),
                                   "Labeling time (days)" = days)
                                 )
                                
tiff(file = "wt_fad_heatmap.png", width = 300, height = 500, units = "px",res = 100)
Heatmap(as.matrix(wt_fad_merge[!is.na(wt_fad_merge[,c(5:8,11:14)] %>% rowMeans()),c(5:8,11:14)]),
                         cluster_columns = F,
                         show_row_dend = F,show_column_names = F,show_row_names = F,
                         top_annotation = wt_fad_anno
        )
dev.off()



####### boxplot and scatter plot of half-lives #######

png("wt_fad_halflives_boxplot_50.png", height = 800, width = 600,units = "px",res = 300)
ggplot(data = melt(wt_fad_merge[which(wt_fad_merge$WT <50),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

png("wt_fad_halflives_boxplot_all.png", height = 800, width = 600,units = "px",res = 300)
ggplot(data = melt(wt_fad_merge[,c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

png("wt_fad_halflives_scatter_200.png", height = 1200, width = 1200,units = "px",res = 150)
ggscatterstats(wt_fad_merge[which(wt_fad_merge$`WT` <200
                                  & wt_fad_merge$FAD <200),],
               x = `WT`,
               y= `FAD`,
               bf.message = F,
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               results.subtitle = T,
               smooth.line.args = list(size = 0.5, color = "blue",alpha = 0.5),
               ggplot.component =
                 list(geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5))
)
dev.off()

png("wt_fad_halflives_scatter_50.png", height = 1200, width = 1200,units = "px",res = 150)
ggscatterstats(wt_fad_merge[which(wt_fad_merge$`WT` <50
                                  & wt_fad_merge$FAD <50),],
               x = `WT`,
               y= `FAD`,
               smooth.line.args = list(size = 0.5, color = "red",alpha = 0.5),
               bf.message = F,
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)"),
ggplot.component =
  list(geom_abline(intercept = 0, slope = 1,
                   color = "red", size = 1,
                   alpha = 0.5))
dev.off()

#### compare WT protein half-lives with Nat Com paper #####
nat_com_halflives <- read.delim("nat_com_2018_protein_halflives.txt")
names(nat_com_halflives)[1] <- "UNIPROT1"
nat_com_halflives <- merge(nat_com_halflives,wt_fad_merge[,c(3,18)],by = "UNIPROT1")
names(nat_com_halflives)[5:6] <- c("Nat Com", "This work")

tiff("nat_com_compararison.tiff", height = 1200, width = 1200,units = "px",res = 300)
ggscatterstats(nat_com_halflives,
               x = `Nat Com`,
               y= `This work`,
               smooth.line.args = list(size = 1, color = "blue",
                                       alpha = 0.5,method = lm),
               bf.message = F,
               xlab = "Nat Com Half-lives (days)",
               ylab = "This work Half-lives (days)",
               results.subtitle = T,
               marginal = F,
ggplot.component =
  list(geom_abline(intercept = 0, slope = 1,
                   color = "red", size = 1,
                   alpha = 0.5)))
dev.off()

################### slower and faster subset ##########################
### define slower and faster subset ###
proteome_slower <- wt_fad_merge[which(wt_fad_merge$`log2fc z` > 1 & !is.na(wt_fad_merge$GN)),] 
write_delim(as.data.frame(proteome_slower[,1:2]),"proteome_slower.txt",col_names = T)

png("proteome_slower_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(proteome_slower[which(proteome_slower$WT <200),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

proteome_faster <- wt_fad_merge[which(wt_fad_merge$`log2fc z` < -2 & !is.na(wt_fad_merge$GN)),]
write_delim(as.data.frame(proteome_faster[,1:2]),"proteome_faster.txt",col_names = F)

png("proteome_faster_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(proteome_faster[which(proteome_faster$WT <200),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

### clusterProfiler analysis of slower and faster proteome ###
proteome_slower_GO <- enrichGO(proteome_slower$GN,`org.Mm.eg.db`,keyType = "SYMBOL",
                               ont = "ALL",pvalueCutoff = 0.01)
proteome_slower_GO <- pairwise_termsim(proteome_slower_GO)
write_xlsx(as.data.frame(proteome_slower_GO@result),"proteome_slower_GO_result.xlsx",format_headers = T)
png("proteome_slower_GO1.png", height = 1200, width = 1200,units = "px",res = 150)
dotplot(proteome_slower_GO,showCategory = 15,label_format = 50) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), 
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1)))
dev.off()
png("proteome_slower_GO2.png", height = 1200, width = 1200,units = "px",res = 150)
emapplot(proteome_slower_GO, cex_label_category=.8, cex_line=.5,showCategory = 15) + coord_cartesian() +
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
dev.off()

proteome_faster_GO <- enrichGO(proteome_faster$GN,`org.Mm.eg.db`,keyType = "SYMBOL",
                               ont = "ALL",pvalueCutoff = 0.01)
proteome_faster_GO <- pairwise_termsim(proteome_faster_GO)
write_xlsx(as.data.frame(proteome_faster_GO@result),"proteome_faster_GO_result.xlsx",format_headers = T)
png("proteome_faster_GO1.png", height = 1200, width = 1200,units = "px",res = 150)
dotplot(proteome_faster_GO,showCategory = 15,label_format = 50) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), 
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1)))
dev.off()
png("proteome_faster_GO2.png", height = 1200, width = 1200,units = "px",res = 150)
emapplot(proteome_faster_GO, cex_label_category=.8, cex_line=.5,showCategory = 15) + coord_cartesian() +
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
dev.off()

################### GSEA analysis of ranked dataset ##########################
# ranked list by mean halfives of wt and fad
mean_rank <- wt_fad_merge[order(wt_fad_merge$mean.halflife,decreasing = T),]
mean_rank <- mean_rank[,c(17:21)]
mean_rank$GeneID2 <- as.character(mean_rank$GeneID2)
mean_rank <- mean_rank[complete.cases(mean_rank$GeneID2),]
mean_rank <- distinct(mean_rank,GeneID2,.keep_all = T)
mean_rank2 <- mean_rank$mean.halflife
names(mean_rank2) <-mean_rank$GeneID2
head(mean_rank2)

mean_gseGO <- gseGO(mean_rank2,ont = "ALL",OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",pvalueCutoff = 0.05)
head(mean_gseGO@result)
mean_gseGO2 <- arrange(mean_gseGO, desc(abs(NES)))

png("mean_gseGO.png", height = 1200, width = 800,units = "px",res = 100)
gseaplot2(title = "GSEA Gene Ontology",
          geneSetID = c("GO:0098800","GO:0005746","GO:0005743","GO:0005201"),
          mean_gseGO,1:10, pvalue_table=F, base_size=14)
dev.off()

mean_gseKEGG <- gseKEGG(mean_rank2,organism = "mmu",pvalueCutoff = 0.05)
head(mean_gseKEGG@result)
mean_gseKEGG2 <- arrange(mean_gseKEGG, desc(abs(NES)))
png("mean_gseKEGG.png", height = 1200, width = 800,units = "px",res = 100)
gseaplot2(title = "GSEA KEGG",mean_gseKEGG2,1:10, pvalue_table=F, base_size=14)
dev.off()

mean_gseWP <- gseWP(mean_rank2,organism = "Mus musculus")
head(mean_gseWP@result)
mean_gseWP2 <- arrange(mean_gseWP, desc(abs(NES)))
png("gseWP.png", height = 1200, width = 800,units = "px",res = 100)
gseaplot2(title = "GSEA WikiPathway",mean_gseWP2,1:10, pvalue_table=F, base_size=14)
dev.off()

remove(mean_gseGO2,mean_gseKEGG2,mean_gseWP2)

# ranked list by log2FCz
if(FALSE) {
z_rank <- wt_fad_merge[order(wt_fad_merge$`log2fc z`,decreasing = T),]
z_rank <- z_rank[,c(1,16)]
z_rank$Uniprot <- sapply(strsplit(z_rank$Uniprot, "\\|"),"[[",2)
z_rank$GeneID <- mapIds(org.Mm.eg.db,keys = z_rank$Uniprot,
                           keytype = "UNIPROT",column = "ENTREZID",)
z_rank <- z_rank[complete.cases(z_rank$GeneID),]
z_rank <- distinct(z_rank,GeneID,.keep_all = T)
z_rank2 <- z_rank$`log2fc z`
names(z_rank2) <-z_rank$GeneID
head(z_rank2)
z_gseGO <- gseGO(z_rank2,ont = "ALL",OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",pvalueCutoff = 0.1)
z_gseKEGG <- gseKEGG(z_rank2,organism = "mmu",pvalueCutoff = 0.1)
z_gseWiki <- gseWP(z_rank2,organism = "Mus musculus")
}



####### compare proteins in subproteomes ##########

# protein protein target in Neuron paper
Neuron_target <- c("App","Slit2","Sfrp1","Smoc1","Htra1","Mdk", "Ntn1",
                   "Cthrc1","Ntn3","Slt3","Lsp1","C4b","Clu","Olfml3","Icam1")

Neuron_target_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% c(Neuron_target)),]

png("Neuron_target_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(Neuron_target_halflives[which(Neuron_target_halflives$WT <200),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

png("Neuron_target_scatter.png", height = 1200, width = 1200,units = "px",res = 300)
ggscatterstats(Neuron_target_halflives,
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               #label.expression = `log2fc of halflife (FAD vs WT)` > 0.5 
               #|`log2fc of halflife (FAD vs WT)` < -0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "white",alpha = 0),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(xlim(c(0,12.5)),
                      ylim(c(0,12.5)),
                      geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))
dev.off()



# protein up but mRNA no change in FAD v.s. WT
protein_top_up <- read.delim("../ProTransAnalysis/top_up_inconsistent_genes.txt",header = F)

protein_top_up_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% c(protein_top_up$V1)),]

tiff("protein_top_up_boxplot.tiff", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(protein_top_up_halflives[which(protein_top_up_halflives$WT <200),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

tiff("protein_top_up_scatter.tiff", height = 1200, width = 1200,units = "px",res = 200)
ggscatterstats(protein_top_up_halflives,
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               label.expression = `log2fc of halflife (FAD vs WT)` > 0,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0.5, color = "blue",
                                       alpha = 0.5, method = lm),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(xlim(c(0,15)),geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))
dev.off()

ggplot(protein_top_up_halflives) +
  geom_point(aes(x = `WT`,y= `FAD`)) +
  xlab("WT Half-lives (days)") +
  ylab("5XFAD Half-lives (days)") +
  xlim(c(0,15)) + 
  geom_abline(intercept = 0, slope = 1,
              color = "red", size = 1,
              alpha = 0.5)



# plaque proteome subset from 8m cutoff
plaque_proteome_8m <- read.delim("../plaque proteome/plaque_8m_12m_top.txt",header = F)
plaque_proteome_8m <- plaque_proteome_8m[-1,]
plaque_proteome_8m_halflives <- wt_fad_merge[which(wt_fad_merge$GN 
                                                       %in% 
                                                         c(plaque_proteome_8m)),]
tiff("plauqe_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(melt(plaque_proteome_8m_halflives[,c(3,9)] %>% 
              filter(WT < 200,FAD < 200),
            value.name = "Half-life (days)", 
            variable.name = "mice type")) +
  geom_boxplot(aes(x = `mice type`,y= `Half-life (days)`,fill= `mice type`))
dev.off()

tiff("plaque_scatter.png", height = 1200, width = 1200,units = "px",res = 200)
ggscatterstats(plaque_proteome_8m_halflives%>% 
                 filter(WT < 200,FAD < 200),
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               label.expression = `log2fc z` >0.4,
               results.subtitle = F,
               marginal = F,
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               smooth.line.args = list(size = 1,
                                       color = "blue",
                                       alpha = 0.5,
                                       method = lm),
               ggplot.component =list(geom_abline(intercept = 0,
                                                  slope = 1,
                                                  color = "red",
                                                  size = 1,
                                                  alpha = 0.5)))
dev.off()

# overlap of protein-up conor and plaque


overlap_up_halflives <- wt_fad_merge[which((wt_fad_merge$GN %in% c(protein_up$V1))
                                            &(wt_fad_merge$GN %in% c(plaque_proteome_8m))),]
write_delim(as.data.frame(overlap_up_halflives[,2]),"overlap.txt",col_names = T)

tiff("overlap_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(melt(overlap_up_halflives[,c(3,9)] %>% 
              filter(WT < 200,FAD < 200),
            value.name = "Half-life (days)", 
            variable.name = "mice type")) +
  geom_boxplot(aes(x = `mice type`,y= `Half-life (days)`,fill= `mice type`))
dev.off()

tiff("overlap_scatter.png", height = 1200, width = 1200,units = "px",res = 200)
ggscatterstats(overlap_up_halflives%>% 
                 filter(WT < 200,FAD < 200),
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               label.expression = `log2fc z` >0.4,
               results.subtitle = F,
               marginal = F,
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               smooth.line.args = list(size = 1,
                                       color = "blue",
                                       alpha = 0.5,
                                       method = lm),
               ggplot.component =list(geom_abline(intercept = 0,
                                                  slope = 1,
                                                  color = "red",
                                                  size = 1,
                                                  alpha = 0.5)))
dev.off()

png("overlap_heatmap.png", height = 800, width = 1200,units = "px",res = 125)
Heatmap(as.matrix(overlap_up_halflives[,c(5:8,11:14)]),
        cluster_columns = F,cluster_rows = T,
        column_title = "Protein Turnover by Time",
        row_labels = overlap_up_halflives$GN,
        column_names_side = "bottom",column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "L%")) +
  Heatmap(as.matrix(overlap_up_halflives[,c(3,9)]),
          cluster_columns = F,cluster_rows = T,
          column_title = "Protein Half-lives",
          row_labels = overlap_up_halflives$GN,
          column_names_side = "bottom",column_names_rot = 0,
          heatmap_legend_param = list(title = "Half-life (days)"))
dev.off()


# presnapse
presynapse <- c("Stx1b","Vamp1","Vamp2", "Snap25", "Syt1","Snap91", "Bin1","Amph"
                ,"Pip5k1c","Snca","Sh3ggl1", "Syn1","Calm1")

presynapse_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% presynapse),]

tiff("presynapse_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(presynapse_halflives[,c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

tiff("presynapse_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(presynapse_halflives,
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               #label.expression = `log2fc of halflife (5xFAD vs WT)` >0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",alpha = 0),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))                                              
dev.off()

# M42_matrisome in Nature Neuronscience paper
M42_matrisome <- c("Smoc1", "Ntn3","Mdk","Ntn1","App","Cthrc1","SFRP1","Spon1",
                   "Olfml3", "Flt1", "Slit2", "Apoe","Slit1", "Dag1","Renbp",
                   "Tmeff2", "Bdh2", "Nxph1", 'Ptn', 'Spock1', "Sdc4", "Gpnmb",
                   "Qprt", "Lrp1", "Col25a1","Col11a1", "Gpc5", "Spock2", 'Htra1',
                   "Frzb","Spock3","Ece1")

M42_matrisome_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% M42_matrisome),]

tiff("M42_matrisome_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(M42_matrisome_halflives[,c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

tiff("M42_matrisome_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(M42_matrisome_halflives,
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               #label.expression = `log2fc of halflife (FAD vs WT)` >0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",alpha = 0),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))                                              
dev.off()


# amloid-beta binding

ABB <- read.delim("../databases/GOs/GO0001540 amyloid-beta binding.txt",header = F)

ABB_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% c(ABB$V2)),]

tiff("ABB_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(ABB_halflives[,c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

tiff("ABB_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(ABB_halflives,
               x = `WT`,
               y= `FAD`,
               label.var = GN,
               label.expression = `log2fc of halflife (FAD vs WT)` >0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",alpha = 0),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))                                              
dev.off()
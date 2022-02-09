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
library(fauxnaif)
library(EnhancedVolcano)
library(corrplot)
library(matrixStats)
library(cluster)
library(VennDiagram)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(stringr)

# load mart to species conversion
human_maRt <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_maRt <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# logic function
'%!in%' <- function(x,y)!('%in%'(x,y))

# function to summarize SD and SE for reshaped data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# fucntion to calculate column-wise and row-wise z
col_z_score <- function(x) {
  a <- t((x-colMeans(x)[col(x)]))/colSds(as.matrix(x))
  return(as.data.frame(t(a)))
}

row_z_score <- function(x) {
  a <- t((x-rowMeans(x)))/rowSds(as.matrix(x))[col(x)]
  return(as.data.frame(t(a)))
}

########## dataset processing and overview of turnover profile ##############
wt_fad_turnover <- read.csv("uni_protein_turnover_final_norm_app.csv") %>%
  as.data.frame()
names(wt_fad_turnover)


wt_fad_turnover <- wt_fad_turnover[,-c(1,grep("sig126",names(wt_fad_turnover)),19:22,39)]

names(wt_fad_turnover)[2:31] <- c(paste(rep("wt",10),c("day-0", "day-4", "day-8", "day-16", "day-32"),1:10, sep = "_"),
                                  paste(rep("fad",5),c("day-0", "day-4", "day-8", "day-16", "day-32"),1:5, sep = "_"),
                                  paste(rep("wt",5),c("day-0", "day-4", "day-8", "day-16", "day-32"),11:15, sep = "_"),
                                  paste(rep("fad",10),c("day-0", "day-4", "day-8", "day-16", "day-32"),6:15, sep = "_"))
names(wt_fad_turnover)


# processing of app anova analysis
App <- wt_fad_turnover[which(wt_fad_turnover$Protein.Accession == "Human APP WDSDPSGTK" |
                        wt_fad_turnover$Protein.Accession == "A-beta LVFFAEDVGSNK"|
                        wt_fad_turnover$Protein.Accession == "Mouse App TEEISEVK"),]
App$Protein.Accession <- c("App","A-beta","App")
# generate JUMPt input turnover cut-off
cutoff_matrix2 <- cbind(!is.na(wt_fad_turnover[,32:36]),
                        !is.na(wt_fad_turnover[,32:36]),
                        !is.na(wt_fad_turnover[,37:41]),
                        !is.na(wt_fad_turnover[,32:36]),
                        !is.na(wt_fad_turnover[,37:41]),
                        !is.na(wt_fad_turnover[,37:41]))

App[,2:31] <- App[,2:31]*cutoff_matrix2[1:3,]
App[,2:31][App[,2:31] == 0] <- NA
App_melt <- pivot_longer(App,
                     cols = 2:31,
                     names_to = "sample",
                     values_to = "L%",
                     values_drop_na = T) %>% as.data.frame()
App_melt <- App_melt[which(App_melt$`L%` != 1),]
App_melt$genotype <- factor(sapply(str_split(App_melt$sample,"_"),"[[",1),
                        ordered = T, level = c("wt","fad"))
App_melt$time <- factor(sapply(str_split(App_melt$sample,"_"),"[[",2),
                    ordered = T,
                    level = c("day-0", "day-4", "day-8", "day-16", "day-32"))

App_res.aov2 <- aov(`L%` ~ Protein.Accession*time, data = App_melt)
summary(App_res.aov2)[[1]][1:2,5][1]
summary(App_res.aov2)[[1]][1:2,5][2]


# same-time point cut-off

cutoff_matrix_b1 <- !is.na(wt_fad_turnover[,32:36])
cutoff_matrix_b2 <- !is.na(wt_fad_turnover[,37:41])

cutoff_matrix <- cutoff_matrix_b1*cutoff_matrix_b2
cutoff_matrix <- cbind(cutoff_matrix,cutoff_matrix,cutoff_matrix,
                       cutoff_matrix,cutoff_matrix,cutoff_matrix)

wt_fad_turnover <- wt_fad_turnover[,1:31]
wt_fad_turnover[-1:-3,2:31] <- wt_fad_turnover[-1:-3,2:31]*cutoff_matrix[-1:-3,]

wt_fad_turnover[,2:31][wt_fad_turnover[,2:31] == 0] <- NA



### case study of anova###
if(FALSE) {
wt_fad_turnover_melt <- pivot_longer(wt_fad_turnover[3763,1:31],
                                     cols = 2:31,
                                     names_to = "sample",
                                     values_to = "L%",
                                     values_drop_na = T) %>% as.data.frame()
wt_fad_turnover_melt <- wt_fad_turnover_melt[which(wt_fad_turnover_melt$`L%` != 1),]

wt_fad_turnover_melt$genotype <- factor(sapply(str_split(wt_fad_turnover_melt$sample,"_"),"[[",1),
                                        ordered = T, level = c("wt","fad"))
wt_fad_turnover_melt$time <- factor(sapply(str_split(wt_fad_turnover_melt$sample,"_"),"[[",2),
                                    ordered = T, level = c("day-0","day-4", "day-8", "day-16", "day-32"))
ggplot(wt_fad_turnover_melt) +
  geom_boxplot(aes(x = time, y = `L%`, fill = genotype)) 
wt_fad_turnover_melt <- wt_fad_turnover_melt[which(wt_fad_turnover_melt$time != "day-32"),]

wt_fad_turnover_melt_aov2 <- aov(`L%` ~ genotype*time, data = wt_fad_turnover_melt)
summary(wt_fad_turnover_melt_aov2)
}

### function for two-way ANOVA  ###
two_way_anova <- function(a){
for ( i in 1:nrow(a)) {
melt <- pivot_longer(a[i,],
                     cols = 2:31,
                     names_to = "sample",
                     values_to = "L%",
                     values_drop_na = T) %>% as.data.frame()
melt <- melt[which(melt$`L%` != 1),]
melt$genotype <- factor(sapply(str_split(melt$sample,"_"),"[[",1),
                        ordered = T, level = c("wt","fad"))
melt$time <- factor(sapply(str_split(melt$sample,"_"),"[[",2),
                    ordered = T,
                    level = c("day-0", "day-4", "day-8", "day-16", "day-32"))
if(nrow(melt) > 6) {
res.aov2 <- aov(`L%` ~ genotype*time, data = melt)
a$`genotype anova p value`[i] <- summary(res.aov2)[[1]][1:2,5][1]
a$`time anova p value`[i] <- summary(res.aov2)[[1]][1:2,5][2]
} else {
  a$`genotype anova p value`[i] <- NA
  a$`time anova p value`[i] <- NA
}
}
  return(a)}
wt_fad_turnover <- two_way_anova(wt_fad_turnover)
wt_fad_turnover$`genotype anova p value`[2] <- summary(App_res.aov2)[[1]][1:2,5][1]
wt_fad_turnover$`time anova p value`[2] <- summary(App_res.aov2)[[1]][1:2,5][2]
wt_fad_turnover$`genotype anova BH FDR` <- p.adjust(p = wt_fad_turnover$`genotype anova p value`,
                                      method = "BH")
sum(is.na(wt_fad_turnover$`genotype anova p value`))
sum(wt_fad_turnover$`genotype anova p value` < 0.05, na.rm = T)
sum(wt_fad_turnover$`genotype anova BH FDR` < 0.05, na.rm = T)

jpeg("turnover_twoway_anova_distribution.jpeg", 
     width = 1000, height = 800, res = 120,units = "px")
ggplot(wt_fad_turnover) +
  geom_density(aes(x = -log10(`genotype anova p value`), y = ..count..,
                   fill = "blue"), alpha = 0.4) +
  geom_density(aes(x = -log10(`genotype anova BH FDR`), y = ..count..,
                   fill = "red"), alpha = 0.4) +
  scale_fill_discrete(labels = c("p value","FDR" )) +
  labs(fill = "Two-way ANOVA", y = "count", x = "-log value") +
  geom_vline(xintercept = 1.3) +
  geom_text(aes(x=2.2, label="p/FDR = 0.05", y=4000), 
            colour="red", text=element_text(size=11)) 
dev.off()

# re-read the turnover data and apply JUMPt input cut-off
wt_fad_turnover2 <- read.csv("uni_protein_turnover_final_norm_app.csv") %>%
  as.data.frame()
wt_fad_turnover2 <- wt_fad_turnover2[,-c(1,grep("sig126",names(wt_fad_turnover2)),19:22,39)]

# JUMPt input turnover cut-off
wt_fad_turnover2 <- wt_fad_turnover2[,1:31]

wt_fad_turnover[,2:31] <- wt_fad_turnover2[,2:31]*cutoff_matrix2

wt_fad_turnover[,2:31][wt_fad_turnover[,2:31] == 0] <- NA

names(wt_fad_turnover)[1] <- "Uniprot"

wt_fad_turnover[which(wt_fad_turnover$Uniprot == "Human APP WDSDPSGTK"),
                c(grep("wt_",names(wt_fad_turnover)),
                  grep("anova", names(wt_fad_turnover)))] <- NA
wt_fad_turnover[which(wt_fad_turnover$Uniprot == "A-beta LVFFAEDVGSNK"),
                c(grep("wt_",names(wt_fad_turnover)))] <- NA
remove(wt_fad_turnover2)
App <- wt_fad_turnover[which(wt_fad_turnover$Uniprot == "Human APP WDSDPSGTK" |
                               wt_fad_turnover$Uniprot == "A-beta LVFFAEDVGSNK"|
                               wt_fad_turnover$Uniprot == "Mouse App TEEISEVK"),]
########## dataset processing and overview of JUMPt result ##############
### reading the file and sort the dataset ###
wt_fad_halflives <- read_xlsx("results_SILAC_TMT_JUMPt_0118_w_cutoff_wo_lys_w_norm.xlsx",sheet = "results") %>% as.data.frame()
App_pep <- wt_fad_halflives %>% filter(geneName == "A-beta LVFFAEDVGSNK" |
                                         geneName == "Human APP WDSDPSGTK")
names(App_pep)[c(1:2,8)] <- c("Uniprot","GN","Half-life (days)")
wt_fad_halflives <- wt_fad_halflives %>% filter(geneName != "A-beta LVFFAEDVGSNK",
                                                  geneName != "Human APP WDSDPSGTK")
wt_fad_halflives$Uniprot <- sapply(strsplit(wt_fad_halflives$Var1,'\\$'),"[[",1)
wt_fad_halflives$mice_type <- sapply(strsplit(wt_fad_halflives$Var1,'\\$'),"[[",2)
wt_fad_halflives$mice_type <- factor(wt_fad_halflives$mice_type, levels = c("WT","FAD"),ordered = TRUE)
wt_fad_halflives$GN <- sapply(strsplit(wt_fad_halflives$geneName,'\\$'),"[[",1)
names(wt_fad_halflives)[c(8,12)] <- c("Half-life (days)","mice type")
wt_fad_halflives <- wt_fad_halflives[-c(grep("HUMAN",wt_fad_halflives$Var1),
                                        grep("cu\\|",wt_fad_halflives$Var1),
                                        grep("co\\|",wt_fad_halflives$Var1)),]
wt_fad_halflives <- wt_fad_halflives %>% mutate(relative_residual_error = residual_error/`Half-life (days)`)
shortlived_halflife <- knots(ecdf(wt_fad_halflives[which(wt_fad_halflives$`Half-life (days)`<1),
                      "Half-life (days)"]))
plot(ecdf(wt_fad_halflives[which(wt_fad_halflives$`Half-life (days)`<1),
                           "Half-life (days)"]))
min_halflife <- shortlived_halflife[round(length(shortlived_halflife)*0.994,0)]
print(min_halflife)
wt_fad_halflives$`Half-life (days)`[wt_fad_halflives$`Half-life (days)` <min_halflife] <- min_halflife

### manage data to seperate and then merge WT amd FAD datasets ###
wt_halflives <- wt_fad_halflives[which(wt_fad_halflives$`mice type` == "WT"),]
wt_halflives <- wt_halflives[,c("Uniprot","GN","time_point_1","time_point_2",
                                "time_point_3","time_point_4","time_point_5",
                                "Half-life (days)")]
wt_halflives <- rbind(wt_halflives,App_pep[,1:8])
wt_halflives[which(wt_halflives$GN == "Human APP WDSDPSGTK"),3:8] <- NA
wt_halflives[which(wt_halflives$GN == "A-beta LVFFAEDVGSNK"),3:8] <- NA
fad_halflives <- wt_fad_halflives[which(wt_fad_halflives$`mice type` == "FAD"),]
fad_halflives <- fad_halflives[,c("Uniprot","GN","time_point_1","time_point_2",
                                  "time_point_3","time_point_4","time_point_5",
                                  "Half-life (days)")]
fad_halflives <- rbind(fad_halflives,App_pep[,1:8])


wt_fad_merge <- merge(wt_halflives,fad_halflives,by = c("Uniprot","GN"))
names(wt_fad_merge)[3:7] <- paste("wt",
                                  c("day-0", "day-4", "day-8", "day-16", "day-32"),
                                  "average", sep = "_")
names(wt_fad_merge)[9:13] <- paste("fad",
                                   c("day-0", "day-4", "day-8", "day-16", "day-32"),
                                   "average", sep = "_")
names(wt_fad_merge)[c(8,14)] <- c("WT Half-life (days)", "5xFAD Half-life (days)")

remove(wt_halflives,fad_halflives)
wt_fad_merge <- merge(wt_fad_turnover,wt_fad_merge,by = c("Uniprot"))

wt_fad_merge[1,"WT Half-life (days)"] = 
  mean(c(wt_fad_merge[2:3,"WT Half-life (days)"],
         wt_fad_merge[2:3,"5xFAD Half-life (days)"]),
       na.rm = T)
### calculate the log2 fold change of half-lives between FAD and WT mice ###
wt_fad_merge <- mutate(wt_fad_merge,
                       `day-4_ratio` = log2(`fad_day-4_average`/`wt_day-4_average`))
wt_fad_merge <- mutate(wt_fad_merge,
                       `day-8_ratio` = log2(`fad_day-8_average`/`wt_day-8_average`))
wt_fad_merge <- mutate(wt_fad_merge,
                       `day-16_ratio` = log2(`fad_day-16_average`/`wt_day-16_average`))
wt_fad_merge <- mutate(wt_fad_merge,
                       `day-32_ratio` = log2(`fad_day-32_average`/`wt_day-32_average`))

wt_fad_merge <- mutate(wt_fad_merge,
                       `Delta half-life (FAD-WT)` = `5xFAD Half-life (days)`-`WT Half-life (days)`
                       )
ggplot(wt_fad_merge) +
  geom_density(aes(x= `Delta half-life (FAD-WT)`, y = ..count..)) +
  xlim(c(-20,20))

wt_fad_merge <- mutate(wt_fad_merge,
                       `Half-life log2FC (FAD/WT)` = log2(`5xFAD Half-life (days)`/`WT Half-life (days)`)
                       )
ggplot(wt_fad_merge) +
  geom_density(aes(x= `Half-life log2FC (FAD/WT)`, y = ..count..)) +
  xlim(c(-5,5))

wt_fad_merge$`Relative half-life change (%)` <- NA
for(i in c(1,3:nrow(wt_fad_merge))){
    if(wt_fad_merge$`WT Half-life (days)`[i] > wt_fad_merge$`5xFAD Half-life (days)`[i]){
    wt_fad_merge$`Relative half-life change (%)`[i] <- 100*wt_fad_merge$`Delta half-life (FAD-WT)`[i]/wt_fad_merge$`WT Half-life (days)`[i]
} else {
  wt_fad_merge$`Relative half-life change (%)`[i] <- 100*wt_fad_merge$`Delta half-life (FAD-WT)`[i]/wt_fad_merge$`5xFAD Half-life (days)`[i]
} 
}

ggplot(wt_fad_merge) +
  geom_density(aes(x= `Relative half-life change (%)`, y = ..count..)) +
  xlim(c(-100,100))

wt_fad_merge <- mutate(wt_fad_merge,
                       mean_halflife = (`5xFAD Half-life (days)`+`WT Half-life (days)`)/2)
wt_fad_merge <- mutate(wt_fad_merge,
                       available_wt_point = 4-is.na(wt_fad_merge[,c("wt_day-4_average",
                                                                  "wt_day-8_average",
                                                                  "wt_day-16_average",
                                                                  "wt_day-32_average")]) %>% 
                         rowMeans()*4)
wt_fad_merge <- mutate(wt_fad_merge,
                       available_fad_point = 4-is.na(wt_fad_merge[,c("fad_day-4_average",
                                                                  "fad_day-8_average",
                                                                  "fad_day-16_average",
                                                                  "fad_day-32_average")]) %>% 
                         rowMeans()*4)

wt_fad_merge <- mutate(wt_fad_merge,
                       available_both_point = (1*!is.na(wt_fad_merge[,c("day-4_ratio",
                                                                      "day-8_ratio",
                                                                      "day-16_ratio",
                                                                      "day-32_ratio")])) %>% rowMeans()*4)
names(wt_fad_merge)



### Gene ID mapping ##
wt_fad_merge$UNIPROT1 <- c("A-beta LVFFAEDVGSNK",
                           "Human APP WDSDPSGTK",
                           "Mouse App TEEISEVK",
                           sapply(strsplit(wt_fad_merge$Uniprot[-1:-3], "\\|"),"[[",2))
wt_fad_merge$UNIPROT2 <- as.character(mapIds(org.Mm.eg.db,keys = wt_fad_merge$GN,
                                keytype = "SYMBOL",column = "UNIPROT"))
wt_fad_merge$GeneID1 <- mapIds(org.Mm.eg.db,keys = wt_fad_merge$UNIPROT1,
                           keytype = "UNIPROT",column = "ENTREZID")
wt_fad_merge$GeneID2 <- mapIds(org.Mm.eg.db,keys = wt_fad_merge$UNIPROT2,
                               keytype = "UNIPROT",column = "ENTREZID")

write_xlsx(wt_fad_merge, "wt_fad_merge.xlsx")

########### overview of half-life dataset ############
### boxplot to overview both half-lives of WT and FAD ###
wt_fad_halflives2 <- wt_fad_halflives
wt_fad_halflives2$`Half-life (days)`[wt_fad_halflives2$`Half-life (days)` >35] <- 36

jpeg("wt_fad_halflives_histogram.jpeg",
     width = 1200, height = 500, units = "px", res = 150)
ggplot(wt_fad_halflives2) +
  geom_histogram(aes(x= `Half-life (days)`,fill = `mice type`),
                 alpha = 0.5,size = 0.01, color = "black",
                 breaks = c(seq(0,35,by=1),36),
                 position = "dodge", show.legend = F)+
  scale_fill_manual(values = rainbow(2)) +
  scale_x_continuous(limits=c(0, 36),
                     breaks=c(seq(0, 35, by=1),36),
                     labels=c(seq(0,35, by=1),"35+"))+
  scale_y_continuous(breaks=c(seq(0, 2000, by=100)))+
  labs(fill = "mice type", y = "Count", x = "Half-life (days)") +
  theme_classic()
dev.off()

jpeg("wt_fad_halflives_density_plot.jpeg",
     width = 2000, height = 1200, units = "px", res = 200)
ggplot(wt_fad_halflives2) +
  geom_density(aes(x= `Half-life (days)`,fill = `mice type`, y = ..count..),
                 alpha = 0.5,size = 0.01, color = "black",
                 breaks = c(seq(0,35,by=1),36),
                 position = "dodge", show.legend = F)+
  scale_fill_manual(values = rainbow(2)) +
  scale_x_continuous(limits=c(0, 36),
                     breaks=c(seq(0, 35, by=1),36),
                     labels=c(seq(0,35, by=1),"35+"))+
  scale_y_continuous(breaks=c(seq(0, 2000, by=100)))+
  labs(fill = "mice type", y = "Count", x = "Half-life (days)") +
  theme_classic()
dev.off()



######### data QC ##############
jpeg("cummulative plot of half-lives.jpeg",
     width = 500, height = 1200, units = "px", res = 150)
par(mfrow = c(3, 1))
plot(ecdf(wt_fad_halflives[,
                           "Half-life (days)"]),
     main = "Cumulative plot (all)",
     xlab = "Half-life (days)") 
plot(ecdf(wt_fad_halflives[which(wt_fad_halflives$`Half-life (days)`<10),
                           "Half-life (days)"]),
     main = "Cumulative plot (Half-life < 10 days)",
     xlab = "Half-life (days)") 
plot(ecdf(wt_fad_halflives[which(wt_fad_halflives$`Half-life (days)`<1),
                           "Half-life (days)"]),
     main = "Cumulative plot (Half-life < 1 days)",
     xlab = "Half-life (days)") 
dev.off()

jpeg("scatter plot of residual errors.jpeg",
     width = 500, height = 1200, units = "px", res = 150)
par(mfrow = c(2,1))
plot(x= wt_fad_halflives$`Half-life (days)`, y = wt_fad_halflives$`residual_error`)
plot(x= wt_fad_halflives$`Half-life (days)`, y = wt_fad_halflives$`relative_residual_error`)
dev.off()

jpeg("scatter plot of residual errors samll range.jpeg",
     width = 500, height = 1200, units = "px", res = 150)
par(mfrow = c(2,1))
plot(x= wt_fad_halflives$`Half-life (days)`, y = wt_fad_halflives$`residual_error`,
     xlim = c(0,1))
plot(x= wt_fad_halflives$`Half-life (days)`, y = wt_fad_halflives$`relative_residual_error`,
     xlim = c(0,1))
dev.off()


turnover_rario <- melt(wt_fad_merge[,c("day-4_ratio",
                                       "day-8_ratio",
                                       "day-16_ratio",
                                       "day-32_ratio")],
                       variable.name = "Days",
                       value.name = "Log2FC(light%, FAD/WT)")


jpeg("turnover_rario.jpeg", 
     width = 1200, height = 800, res = 150,units = "px")
ggplot(turnover_rario) +
  geom_density(aes(x= `Log2FC(light%, FAD/WT)`, y = ..scaled.., color = `Days`),
               fill = "grey", alpha = 0.2)
dev.off()

half_life_diff <- melt(wt_fad_merge[,c("Delta half-life (FAD-WT)",
                                       "Half-life log2FC (FAD/WT)",
                                       "Relative half-life change (%)")],
                       variable.name = "Method",
                       value.name = "Half-life difference")


jpeg("half_life_diff_all.jpeg", 
     width = 1000, height = 800, res = 120,units = "px")
ggplot(half_life_diff) +
  geom_density(aes(x = `Half-life difference`, y = ..scaled..,fill = Method),
               size = 1, alpha = 0.5) +
  #xlim(c(-10,10)) +
  xlab("Half-life difference") +
  ylab("Scaled density") 
dev.off()

jpeg("half_life_diff_zoomin.jpeg", 
     width = 1000, height = 800, res = 120,units = "px")
ggplot(half_life_diff) +
  geom_density(aes(x = `Half-life difference`, y = ..scaled..,fill = Method),
               size = 1, alpha = 0.5) +
  xlim(c(-100,100)) +
  xlab("Half-life difference") +
  ylab("Scaled density")
dev.off()





###### overview of L% #
if(FALSE) {
# selected L% heatmap of both mice types
days = colorRamp2(c(4, 32), c("white", "blue"))
wt_fad_anno <- HeatmapAnnotation(Phenotype = factor(c("WT","WT","WT","WT","5XFAD","5XFAD","5XFAD","5XFAD"),
                                                    ordered = TRUE,levels = c("WT","5XFAD")),
                                 "Labeling time (days)" = c(4,8,16,32,4,8,16,32),
                                 col = list(Phenotype = setNames(rainbow(2),c("WT","5XFAD")),
                                   "Labeling time (days)" = days)
                                 )
                                
jpeg(file = "wt_fad_heatmap.png", width = 300, height = 500, units = "px",res = 100)
Heatmap(as.matrix(wt_fad_merge[!is.na(wt_fad_merge[,c(5:8,11:14)] %>% rowMeans()),c(5:8,11:14)]),
                         cluster_columns = F,
                         show_row_dend = F,show_column_names = F,show_row_names = F,
                         top_annotation = wt_fad_anno
        )
dev.off()
}

### compare WT protein half-lives with Nat Com paper #####
nat_com_halflives <- read.delim("nat_com_2018_protein_halflives.txt")
names(nat_com_halflives)[1] <- "UNIPROT1"
nat_com_halflives_merge <- merge(nat_com_halflives,
                           wt_fad_merge[,c("WT Half-life (days)",
                                           "UNIPROT1")],
                           by = "UNIPROT1")
names(nat_com_halflives_merge)[19] <- c("This work")

nat_com_halflives_cor <- cor(nat_com_halflives_merge[,5:19],use = "na.or.complete")

jpeg("nat_com_corr.jpeg", height = 1200, width = 1200,units = "px",res = 300)
corrplot(nat_com_halflives_cor,method = 'color',diag = F,
         is.corr =F,tl.pos =F,type = 'lower')
dev.off()



nat_com_halflives_merge$brain_average <- 
  rowMeans(nat_com_halflives_merge[,5:13],
           na.rm = T)
nat_com_halflives_merge$cortex_average <- 
  rowMeans(nat_com_halflives_merge[,5:6],
           na.rm = T)

jpeg("nat_com_compararison.jpeg", height = 1200, width = 1200,units = "px",res = 300)
ggscatterstats(nat_com_halflives_merge,
               x = `Protein.half.life.in.brain.cortex.homogenate.control..days.`,
               y= `This work`,
               smooth.line.args = list(size = 1, color = "blue",
                                       alpha = 0.5,method = lm,
                                       formula = y ~ x+0),
               bf.message = F,
               xlab = "Half-lives (days) in Fornasiero's work",
               ylab = "Half-lives (days) in this work ",
               results.subtitle = T,
               marginal = F,
ggplot.component =
  list(geom_abline(intercept = 0, slope = 1,
                   color = "red", size = 1,
                   alpha = 0.5),
       theme_classic()))
dev.off()

venn.diagram(
  x = list(wt_fad_merge$UNIPROT1,
           nat_com_halflives$UNIPROT1[
             which(!is.na(nat_com_halflives$Protein.half.life.in.brain.cortex.homogenate.control..days.))]),
  category.names = c("This work" , "Fornasiero's work" ),
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
  rotation.degree = 45,
)

### compare wt half-lives with brain concentration estimated by APEX

brain_protein_concentration <- read_excel("Brain_protein_concentrations.xlsx", skip = 2)
names(brain_protein_concentration)[2] <- "Uniprot"
brain_protein_concentration <- merge(brain_protein_concentration,
                                     wt_fad_merge[,c("Uniprot",
                                                     "WT Half-life (days)")],
                                     by = "Uniprot")
brain_protein_concentration <- brain_protein_concentration[,-7]
brain_conc_cor <- cor(brain_protein_concentration[,5:7],
                      use = "na.or.complete")

jpeg("brain_conc_halflife_corlot.jpeg",
     width = 900, height = 900, units = "px", res = 200)
corrplot(brain_conc_cor,method = 'color',diag = F,
         is.corr =F,tl.pos =F,type = 'lower')
dev.off()


### compare wt half-lives with protein properties

gravy <- read_excel("mouse_proteome_GRAVY.xlsx")
names(gravy)[1] <- "UNIPROT1"
pi_mw <- read.delim("mouse_proteome_pi_mw.txt",header = F)
names(pi_mw)[2:4] <- c("UNIPROT1", "PI","Mw") 

protein_property <- merge(wt_fad_merge[,c("UNIPROT1","WT Half-life (days)"),],
                          gravy[,1:2] %>% as.data.frame(), 
                          by = "UNIPROT1", all = T)
protein_property <- merge(protein_property, pi_mw[,2:4],
                          by = "UNIPROT1", all = T)
protein_property$GRAVY <- as.numeric(protein_property$GRAVY)
protein_property$PI <- as.numeric(protein_property$PI)
protein_property_cor <- cor(protein_property[,2:5],
                            use = "na.or.complete")
jpeg("protein_property_halflife_corlot.jpeg",
     width = 900, height = 900, units = "px", res = 200)
corrplot(protein_property_cor,method = 'color',diag = F,
         is.corr =F,tl.pos =F,type = 'lower')
dev.off()


### compare wt half-lives by cellular location
location <- read_excel("protein cellular location.xlsx")



human2mouse <- getLDS(attributes = "uniprotswissprot", 
                         filters = "uniprotswissprot", 
                         values = location$Uniprot, 
                         mart = human_maRt, 
                         attributesL = "uniprotswissprot", 
                         martL = mouse_maRt)
names(human2mouse) <- c("Uniprot","UNIPROT1")
location <- merge(location[,3:7], human2mouse, by = "Uniprot")
location <- merge(location, by = "UNIPROT1",
                  wt_fad_merge[,c("UNIPROT1","WT Half-life (days)")])
location_halflives <- pivot_longer(location,cols = 3:6,names_to = "location_n",
                                   values_to = "location")
location_halflives <- location_halflives %>% filter(location != "NA")

location_unique <- unique(location_halflives$`location`)

for (i in 1:length(location_unique)) {
  location_halflives[which(location_halflives$location == location_unique[i]),5] <-
  paste(location_unique[i]," (N = ",
        sum(location_halflives$location == location_unique[i]),")",
        sep = "")
}

location_halflives$count <- sapply(strsplit(location_halflives$location,"="),"[[",2)
location_halflives$count <- sapply(strsplit(location_halflives$count,")"),"[[",1) %>%
  as.numeric()

location_halflives <- location_halflives %>% filter(count > 10)

location_halflives$`WT Half-life (days)` <- log2(location_halflives$`WT Half-life (days)`)

jpeg("location_halflives.jpeg",
     width = 1200, height = 900, units = "px", res = 200)
ggplot(location_halflives) +
  geom_boxplot(aes(y = reorder(location,`WT Half-life (days)`, FUN = median),
                   x = `WT Half-life (days)`, fill = location),
               show.legend = F) +
  #xlim(c(0,40)) +
  xlab("Log2-transformed half-life (days)") +
  ylab("Subcellular location")
dev.off()

############ turnover of specific proteins or peptides #######
wt_fad_merge_melt <- pivot_longer(wt_fad_merge[,c(1:31,35)], cols = 2:31, 
                                  names_to = "sample name",
                                  values_to = "turnover")
wt_fad_merge_melt$phenotype <- sapply(strsplit(wt_fad_merge_melt$`sample name`,"_"),"[[",1)
wt_fad_merge_melt$time_1 <- sapply(strsplit(wt_fad_merge_melt$`sample name`,"_"),"[[",2)
wt_fad_merge_melt$time_1 <- factor(wt_fad_merge_melt$time_1, ordered = T,
                                   levels = c("day-0","day-4",
                                              "day-8","day-16","day-32"))
wt_fad_merge_melt$time_2 <- rep(c(0,4,8,16,32),nrow(wt_fad_merge_melt)/5)



App <- data_summary(wt_fad_merge_melt %>%
                      filter(Uniprot %in% c("A-beta LVFFAEDVGSNK",
                                       "Human APP WDSDPSGTK",
                                       "Mouse App TEEISEVK")),
                    varname="turnover",
                    groupnames=c("GN","phenotype","time_2"))
App <- App %>% filter(sd != "NA")

jpeg("App.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(App,
       aes(x=time_2, y=turnover, shape = GN, color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

# slower proteins
Apoe <- data_summary(wt_fad_merge_melt %>%
                                    filter(GN == "Apoe"),
                                  varname="turnover",
                                  groupnames=c("phenotype","time_2"))
jpeg("Apoe.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Apoe,
        aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

Sfrp1 <- data_summary(wt_fad_merge_melt %>%
                       filter(GN == "Sfrp1"),
                     varname="turnover",
                     groupnames=c("phenotype","time_2"))
jpeg("Sfrp1.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Sfrp1,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

Slit2 <- data_summary(wt_fad_merge_melt %>%
                        filter(GN == "Slit2"),
                      varname="turnover",
                      groupnames=c("phenotype","time_2"))
jpeg("Slit2.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Slit2,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

Itm2c <- data_summary(wt_fad_merge_melt %>%
                        filter(GN == "Itm2c"),
                      varname="turnover",
                      groupnames=c("phenotype","time_2"))
jpeg("Itm2c.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Itm2c,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

# faster proteins
Tulp2 <- data_summary(wt_fad_merge_melt %>%
                        filter(GN == "Tulp2"),
                      varname="turnover",
                      groupnames=c("phenotype","time_2"))
jpeg("Tulp2.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Tulp2,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

Ctnnd1 <- data_summary(wt_fad_merge_melt %>%
                        filter(Uniprot == "sp|P30999|CTND1_MOUSE"),
                      varname="turnover",
                      groupnames=c("phenotype","time_2"))
jpeg("Ctnnd1.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Ctnnd1,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()


Vwa1 <- data_summary(wt_fad_merge_melt %>%
                        filter(GN == "Vwa1"),
                      varname="turnover",
                      groupnames=c("phenotype","time_2"))
jpeg("Vwa1.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Vwa1,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()

Agrn <- data_summary(wt_fad_merge_melt %>%
                         filter(Uniprot == "sp|A2ASQ1|AGRIN_MOUSE"),
                       varname="turnover",
                       groupnames=c("phenotype","time_2"))
jpeg("Agrn.jpeg",
     width = 500, height = 500, units = "px", res = 150)
ggplot(Agrn,
       aes(x=time_2, y=turnover,  color = phenotype)) + 
  geom_line(size = 1,show.legend = F) +
  geom_point(size = 3,show.legend = F) +
  geom_errorbar(aes(ymin=turnover-sd, ymax=turnover+sd), width= 1,
                position="identity",show.legend = F) +
  theme_classic() +
  xlab("Time (days)") +
  ylab("L/(L+H)")
dev.off()
################### slower and faster subset ##########################
### define slower and faster subset ###
wt_fad_merge_ava <- filter(wt_fad_merge,
                           available_wt_point >0,
                           available_fad_point >0)
wt_fad_merge_ava <- rbind(wt_fad_merge[1,],wt_fad_merge_ava)
mean(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
median(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
sd(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
nrow(wt_fad_merge_ava)

jpeg("RHC_histogram.jpeg", 
     width = 1200, height = 800, res = 200,units = "px")
ggplot(wt_fad_merge_ava) +
  geom_histogram(aes(x= `Relative half-life change (%)`),
                 binwidth = 5, position = "dodge",
                 color = "black", fill = "blue", alpha = 0.5) +
  theme_classic()
dev.off()

wt_fad_merge_ava <- filter(wt_fad_merge,
                           available_both_point >0,
                           !is.na(`genotype anova BH FDR`))
wt_fad_merge_ava <- rbind(wt_fad_merge[1,],wt_fad_merge_ava)
wt_fad_merge_ava <- wt_fad_merge_ava %>% group_by(GN) %>% 
  slice_max(order_by = abs(`Relative half-life change (%)`), n = 1, with_ties = F)

mean(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
median(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
sd(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)
nrow(wt_fad_merge_ava)

rhc_sd <- sd(wt_fad_merge_ava$`Relative half-life change (%)`, na.rm = T)*1
print(rhc_sd)
proteome_slower <- wt_fad_merge_ava[which(wt_fad_merge_ava$`Relative half-life change (%)` > rhc_sd & 
                                        !is.na(wt_fad_merge_ava$GN) &
                                        wt_fad_merge_ava$`genotype anova BH FDR` < 0.05),] 
proteome_slower <- proteome_slower[order(proteome_slower$`Relative half-life change (%)`,decreasing = T),]
print(nrow(proteome_slower))

proteome_faster <- wt_fad_merge_ava[which(wt_fad_merge_ava$`Relative half-life change (%)` < 0-rhc_sd & 
                                        !is.na(wt_fad_merge_ava$GN) &
                                        wt_fad_merge_ava$`genotype anova BH FDR` < 0.05),]
proteome_faster <- proteome_faster[order(proteome_faster$`Relative half-life change (%)`),]
print(nrow(proteome_faster))
write_delim(c(proteome_slower$GN,
              proteome_faster$GN,
              "App") %>% 
              as.data.frame(), 
            "turnover_DE.txt",
            col_names = F)

write_delim(proteome_slower, 
            "turnover_DE_slower.txt",
            col_names = T)


 
 
jpeg("wt_fad_halflife_volcano.jpeg", width = 1600, height = 1800, units = "px", res = 200)
EnhancedVolcano(toptable = wt_fad_merge_ava,
                x = "Relative half-life change (%)", 
                y = "genotype anova BH FDR",
                lab = wt_fad_merge_ava$GN,
                ylim = c(0,17), pointSize = 2,
                selectLab = c(proteome_slower$GN[1:20],proteome_faster$GN[1:10]),
                FCcutoff = rhc_sd, pCutoff = 0.05,
                xlab = "Relative half-life change (%)",ylab = "-Log(ANOVA FDR)",
                gridlines.major = F,gridlines.minor = F,
                drawConnectors = T,
                maxoverlapsConnectors = Inf,
                title = "Half-life volcano plot", subtitle = "5xFAD v.s. WT")
dev.off()

jpeg("wt_fad_halflife_volcano2.jpeg", width = 1600, height = 1800, units = "px", res = 200)
EnhancedVolcano(toptable = wt_fad_merge_ava,
                x = "Relative half-life change (%)", 
                y = "genotype anova BH FDR",
                lab = wt_fad_merge_ava$GN,
                ylim = c(0,4), pointSize = 2,
                selectLab = c(proteome_slower$GN[1:20],proteome_faster$GN[1:10]),
                FCcutoff = rhc_sd, pCutoff = 0.05,
                xlab = "Relative half-life change (%)",ylab = "-Log(ANOVA FDR)",
                gridlines.major = F,gridlines.minor = F,
                drawConnectors = T,
                maxoverlapsConnectors = Inf,
                title = "Half-life volcano plot", subtitle = "5xFAD v.s. WT")
dev.off()

###### jump
turnover_de_jumpn <- read.delim("z:/ResearchHome/ClusterHome/kyang2/JUMPn/turnover_DE/intermediate/PPI/foreground.txt_Publication_table_FDRfiltered.txt")

turnover_de_jumpn <- turnover_de_jumpn %>% 
  filter(GeneSetName != "GeneSetName")


turnover_de_jumpn[,c(3:5,7:9)] <- turnover_de_jumpn[,c(3:5,7:9)] %>% 
  as.matrix() %>% 
  as.numeric()
turnover_de_jumpn$Database <- 
  sapply(strsplit(turnover_de_jumpn$GeneSetName, split = "_"), "[[",1)
turnover_de_jumpn <- turnover_de_jumpn %>% 
  mutate(BH = 0-log10(BH))
names(turnover_de_jumpn)[9] <- "-log(BH FDR)"

turnover_de_jumpn$Pathway <-
  str_replace(turnover_de_jumpn$GeneSetName,
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
  str_replace_all(pattern = "_",replacement = " ")

turnover_de_jumpn$Pathway <- turnover_de_jumpn$Pathway %>% tolower() %>% str_to_title()
turnover_de_jumpn$name_length <- nchar(turnover_de_jumpn$Pathway)

turnover_de_jumpn <- turnover_de_jumpn %>%
  filter(`-log(BH FDR)` > 1.3) %>%
  filter(name_length < 40) %>%
  arrange(Database,FoldEnrichment)

turnover_de_jumpn <- turnover_de_jumpn %>%
  filter(Pathway %!in% c("Isoprenoid Metabolic Process",
                         "Amyloid Fibril Formation",
                         "Sulfur Compound Binding",
                         "Glycosaminoglycan Binding",
                         "Heparan Sulfate Proteoglycan Binding",
                         "Signaling By Erbb4",
                         "Hs Gag Biosynthesis"))

turnover_de_jumpn$Pathway <- factor(turnover_de_jumpn$Pathway,
                                    levels = turnover_de_jumpn$Pathway)

jpeg("turnover_de_jumpn.jpeg",
     width = 1200, height = 1200,
     units = "px", res = 200)
ggplot(turnover_de_jumpn) +
  geom_point(aes(x = FoldEnrichment, y = Pathway,
                 color = Database, size = `-log(BH FDR)`)) +
  xlab("Fold of enrichment") +
  ylab("M2-enriched pathways") +
  theme_bw()
dev.off()



####### clusterProfiler analysis of slower and faster proteome ##############
enrich_GO_result <- function(a, group = "group") {
  a <- enrichGO(a,'org.Mm.eg.db',
                keyType = "SYMBOL",ont="ALL",
                pvalueCutoff=1,pAdjustMethod = "BH")
  a <- a@result
  names(a)[4:10] <- paste(group,names(a)[4:10], sep = "_")
  return(a)
}
proteome_slower_GO_result <- enrich_GO_result(proteome_slower$GN,
                                       group = "slower")

proteome_faster_GO_result <- enrich_GO_result(proteome_faster$GN,
                                       group = "faster")

C1M2_GO_result <- enrich_GO_result(c("Ifit3","Atp9a","Apoe",
                                     "Agrn","Gpc1","Slit2"),
                                              group = "C1M2")

proteome_slower_GO_result_select <- filter(proteome_slower_GO_result,
                                           Description == "lysosome")[,"slower_geneID"] %>%
  strsplit('/')
proteome_slower_GO_result_select <- proteome_slower_GO_result_select[[1]]

proteome_slower_faster_GO_result <- merge(proteome_slower_GO_result,
                                          proteome_faster_GO_result,
                                          by= c("ID","Description","ONTOLOGY"), all = T)
proteome_slower_faster_GO_result <- proteome_slower_faster_GO_result[order(proteome_slower_faster_GO_result$slower_p.adjust),
                                                   c(1:3,
                                                     grep("p.adjust",names(proteome_slower_faster_GO_result)),
                                                     grep("geneID",names(proteome_slower_faster_GO_result)))]
proteome_slower_faster_GO_result[4:5] <- log10(proteome_slower_faster_GO_result[4:5])*(-1)
proteome_slower_faster_GO_result[is.na(proteome_slower_faster_GO_result)] <- 0

proteome_slower_faster_GO_result %>% 
  filter(ONTOLOGY == "BP") %>%
  slice_max(order_by = faster_p.adjust,n = 20) 

proteome_slower_faster_GO_top <- c("lysosome",
                                   "phagocytic vesicle",
                                   "endocytic vesicle",
                                   "collagen-containing extracellular matrix",
                                   "inhibitory synapse",
                                   "negative regulation of cell projection organization",
                                   "lipid transport",
                                   "glial cell differentiation",
                                   "glycosaminoglycan binding",
                                   "lipid transporter activity",
                                   "amyloid-beta binding",
                                   "double-stranded RNA binding",
                                   "tau protein binding",
                                   
                                   "extracellular matrix structural constituent",
                                   "growth factor binding",
                                   "synaptic cleft",
                                   "heparan sulfate proteoglycan binding",
                                   "positive regulation of growth")

FDR_GO <- colorRamp2(c(0,2,4),c("white", "steelblue3","steelblue4"))

proteome_slower_faster_GO_result_select <- filter(proteome_slower_faster_GO_result, 
                                         Description %in% proteome_slower_faster_GO_top)

jpeg("slower_faster_GO_heatmap.jpeg", width = 900, height = 900, units = "px", res = 140)
Heatmap(proteome_slower_faster_GO_result_select[,4:5] %>% 
          as.matrix(),
        row_labels = proteome_slower_faster_GO_result_select$Description,
        row_names_side = "right", show_column_names = F,
        cluster_rows = diana(proteome_slower_faster_GO_result_select[,4:5]), 
        cluster_columns = F, show_row_dend = F,
        #row_split = proteome_slower_faster_GO_result_select$ONTOLOGY,
        col = FDR_GO, border = T,row_names_max_width = unit(12, "cm"),
        row_names_gp = gpar(fontsize = 16),
        heatmap_legend_param = list(title = "-LogFDR"))
dev.off()



###### master table of four datasets #########
### mouse master table of proteome,transcritome and turnover ###
pro_trans <- read_xlsx("../ProTransAnalysis/pro_trans.xlsx")
pro_trans <- pro_trans[,c(1:3,7,5:6,13,11:12)]
for(i in 1:nrow(pro_trans)) {
  if(pro_trans$`5XFAD v.s. WT protein log2FC-z`[i] >0){
    pro_trans$`protein change`[i] <- "protein up"
  } else {
    pro_trans$`protein change`[i] <- "protein down"
    }
}
for(i in 1:nrow(pro_trans)) {
  if(pro_trans$`5XFAD v.s. WT mRNA log2FC-z`[i] >0){
    pro_trans$`mRNA change`[i] <- "mRNA up"
  } else {
    pro_trans$`mRNA change`[i] <- "mRNA down"
  }
}

for(i in 1:nrow(pro_trans)){
  if((pro_trans$`5XFAD v.s. WT mRNA log2FC-z`[i])*(pro_trans$`5XFAD v.s. WT protein log2FC-z`[i]) > 0) {
    if(abs(pro_trans$`5XFAD v.s. WT protein log2FC-z`[i] - pro_trans$`5XFAD v.s. WT mRNA log2FC-z`[i]) > 3
    ) {pro_trans$`consistency`[i] <- "inconsistent"
    }else{pro_trans$`consistency`[i] <- "consistent"}
  } else {
    pro_trans$`consistency`[i] <- "inconsistent"
  }}
pro_trans[which(pro_trans$`Protein Accession` == "sp|P05067|A4_HUMAN"),1:3] <- "A-beta LVFFAEDVGSNK"

plaque_proteome <- read_xlsx("../plaque proteome/plaque_proteome.xlsx")
names(plaque_proteome)[1] <- c("Protein Accession")
plaque_proteome <- plaque_proteome[,c(1,4:8)]
for(i in 1:nrow(plaque_proteome)) {
  if(plaque_proteome$`8M and 12M plq v.s. nonplq log2FC-z`[i] >0){
    plaque_proteome$`plaque protein change`[i] <- "protein up"
  } else {
    plaque_proteome$`plaque protein change`[i] <- "protein down"
  }
}

turnover <- wt_fad_merge[,c(1,
                            grep("GN",names(wt_fad_merge)),
                            grep("Half-life",names(wt_fad_merge)),
                            grep("anova",names(wt_fad_merge)),
                            grep("half-life",names(wt_fad_merge)))]
names(turnover)[1] <- "Protein Accession"

turnover_meta_table <- merge(pro_trans,plaque_proteome, by = c("Protein Accession"), all =T)
turnover_meta_table <- merge(turnover_meta_table,turnover,by = c("Protein Accession"), all =T)
turnover_meta_table <- turnover_meta_table[-grep("CON_",turnover_meta_table$`Protein Accession`),]

for(i in 1:nrow(turnover_meta_table)) {
  if(turnover_meta_table$GN.x[i] %>% is.na()){
    if(turnover_meta_table$GN.y[i] %>% is.na()){
      turnover_meta_table$GN.x[i] <- turnover_meta_table$GN[i]
    } else{
    turnover_meta_table$GN.x[i] <- turnover_meta_table$GN.y[i]
    }
  }
}
turnover_meta_table <- turnover_meta_table[,which(names(turnover_meta_table) != "GN.y" & 
                                                    names(turnover_meta_table) != "GN")]
names(turnover_meta_table)[2] <- "GN"

names(turnover_meta_table)[1:3] <- paste("Mouse",names(turnover_meta_table)[1:3], sep = " ")

write_xlsx(turnover_meta_table,"turnover_mata_table_mouse.xlsx")
################## correlation analysis of meta table ########

turnover_meta_cor <- cor(turnover_meta_table[,c("5XFAD v.s. WT protein log2FC-z",
                                                "5XFAD v.s. WT mRNA log2FC-z",
                                                "Relative half-life change (%)"
)],
use = "na.or.complete")
jpeg("protrans_turnover_corlot.jpeg", width = 900, height = 900, units = "px", res = 200)
corrplot(turnover_meta_cor,method = 'color',diag = F,
         is.corr = T,tl.pos =F,type = 'lower')
dev.off()


jpeg("protrans_turnover_boxplot.jpeg", width = 900, height = 900, units = "px", res = 200)
ggplot(turnover_meta_table[!is.na(turnover_meta_table$consistency),]) +
  geom_boxplot(aes(x =`consistency` ,fill =`protein change`,color = `mRNA change`,
                   y = `Relative half-life change (%)`)) +
  #scale_fill_manual(values = rainbow(2)) + 
  ylab("Relative half-life change (%)") +
  ylim(c(-100,100)) +
  geom_hline(yintercept = 0)+
  theme_classic() 
dev.off()


con_pro_up <- filter(turnover_meta_table,
                     consistency == "consistent",
                     `protein change` == "protein up")[,"Relative half-life change (%)"]
con_pro_down <- filter(turnover_meta_table,
                     consistency == "consistent",
                     `protein change` == "protein down")[,"Relative half-life change (%)"]



incon_pro_up <- filter(turnover_meta_table,
                     consistency == "inconsistent",
                     `protein change` == "protein up")[,"Relative half-life change (%)"]
incon_pro_down <- filter(turnover_meta_table,
                       consistency == "inconsistent",
                       `protein change` == "protein down")[,"Relative half-life change (%)"]

incon_pro_up_mrna_up <- filter(turnover_meta_table,
                       consistency == "inconsistent",
                       `protein change` == "protein up",
                       `mRNA change` == "mRNA up")[,"Relative half-life change (%)"]
incon_pro_up_mrna_down <- filter(turnover_meta_table,
                               consistency == "inconsistent",
                               `protein change` == "protein up",
                               `mRNA change` == "mRNA down")[,"Relative half-life change (%)"]
incon_pro_down_mrna_up <- filter(turnover_meta_table,
                         consistency == "inconsistent",
                         `protein change` == "protein down",
                         `mRNA change` == "mRNA up")[,"Relative half-life change (%)"]
incon_pro_down_mrna_down <- filter(turnover_meta_table,
                                 consistency == "inconsistent",
                                 `protein change` == "protein down",
                                 `mRNA change` == "mRNA down")[,"Relative half-life change (%)"]

wilcox.test(x= con_pro_up, y = con_pro_down)[[3]]
wilcox.test(x= incon_pro_up, y = incon_pro_down)[[3]]
wilcox.test(x= incon_pro_up_mrna_up, y = incon_pro_up_mrna_down)[[3]]
wilcox.test(x= incon_pro_down_mrna_up, y = incon_pro_down_mrna_down)[[3]]
wilcox.test(x= incon_pro_up_mrna_down, y = incon_pro_down_mrna_down)[[3]]
wilcox.test(x= incon_pro_up_mrna_up, y = incon_pro_down_mrna_up)[[3]]
wilcox.test(x= incon_pro_up_mrna_up, y = incon_pro_down_mrna_down)[[3]]
wilcox.test(x= incon_pro_up_mrna_down, y = incon_pro_down_mrna_up)[[3]]










######### check the protein half-life by categories defined from protrans data ############
incon_proteins <- read.delim("../ProTransAnalysis/incon_proteins.txt")

turnover_meta_incon <- turnover_meta_table %>% 
  filter(`Mouse Protein Accession` %in% incon_proteins$x)
turnover_meta_incon <- turnover_meta_incon[!is.na(turnover_meta_incon$`Mouse GN`),]
col_z_protein = colorRamp2(c(min(turnover_meta_incon$`5XFAD v.s. WT protein log2FC-z`,na.rm = T)
                             ,0,
                             max(turnover_meta_incon$`5XFAD v.s. WT protein log2FC-z`,na.rm = T)),
                           c("blue", "white", "red"))
col_z_mrna = colorRamp2(c(min(turnover_meta_incon$`5XFAD v.s. WT mRNA log2FC-z`,na.rm = T)
                          ,0,
                          max(turnover_meta_incon$`5XFAD v.s. WT mRNA log2FC-z`,na.rm = T)),
                        c("blue", "white", "red"))
col_z_rhc = colorRamp2(c(-100,0,100),
                        c("blue", "white", "red"))

jpeg("selected_turnover_meta.jpeg", height = 1000, width = 500, units = "px", res = 150)
Heatmap(turnover_meta_incon[,"5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = turnover_meta_incon$`Mouse GN`,
        row_names_side = "left",show_row_dend = F,
        cluster_rows = turnover_meta_incon[,c("5XFAD v.s. WT protein log2FC-z",
                                              "5XFAD v.s. WT mRNA log2FC-z")] %>% diana(),
        show_heatmap_legend = T,
        col = col_z_protein, border = T) +
  Heatmap(turnover_meta_incon[,"5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = T,
          col = col_z_protein, border = T) +
  Heatmap(turnover_meta_incon[,"Relative half-life change (%)"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = T,
          col = col_z_rhc, border = T)
dev.off()

########### check the protein half-life by categories enriched from GO ##############

proteome_slower_faster <- rbind(proteome_slower,proteome_faster)
proteome_slower_faster <- turnover_meta_table %>%
  filter(`Mouse Protein Accession` %in% proteome_slower_faster$Uniprot) 
proteome_slower_faster[1,"Mouse GN"] <- "A-beta"

proteome_slower_faster <- 
  proteome_slower_faster[order(
    proteome_slower_faster$`Relative half-life change (%)`,
    na.last= T,
    decreasing = T),]

col_z_protein = colorRamp2(c(min(proteome_slower_faster$`5XFAD v.s. WT protein log2FC-z`,na.rm = T)
                             ,0,
                             max(proteome_slower_faster$`5XFAD v.s. WT protein log2FC-z`,na.rm = T)),
                           c("blue", "white", "red"))
col_z_mrna = colorRamp2(c(min(proteome_slower_faster$`5XFAD v.s. WT mRNA log2FC-z`,na.rm = T)
                          ,0,
                          max(proteome_slower_faster$`5XFAD v.s. WT mRNA log2FC-z`,na.rm = T)),
                        c("blue", "white", "red"))
col_z_plaque = colorRamp2(c(min(proteome_slower_faster$`8M and 12M plq v.s. nonplq log2FC-z`,na.rm = T)
                          ,0,
                          max(proteome_slower_faster$`8M and 12M plq v.s. nonplq log2FC-z`,na.rm = T)),
                        c("blue", "white", "red"))

jpeg("DE_meta.jpeg", height = 1000, width = 500, units = "px", res = 150)
Heatmap(proteome_slower_faster[,"5XFAD v.s. WT protein log2FC-z"] %>% as.matrix(),
        show_column_names = F,cluster_columns = F,
        row_labels = proteome_slower_faster$`Mouse GN`,
        row_names_side = "left",show_row_dend = F,
        cluster_rows = F,
        show_heatmap_legend = T,heatmap_legend_param = list(title = ""),
        col = col_z_protein, border = T) +
  Heatmap(proteome_slower_faster[,"5XFAD v.s. WT mRNA log2FC-z"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = T,heatmap_legend_param = list(title = ""),
          col = col_z_protein, border = T) +
  Heatmap(proteome_slower_faster[,"Relative half-life change (%)"] %>% as.matrix(),
          show_column_names = F,cluster_columns = F,
          show_row_names = F,cluster_rows = F,
          show_heatmap_legend = T,heatmap_legend_param = list(title = ""),
          col = col_z_rhc, border = T) 
dev.off()
write_xlsx(proteome_slower_faster,"turnover_mata_table_accepted_mouse.xlsx")

proteome_slower_faster_M2 <- proteome_slower_faster %>% 
  filter(`Mouse GN` %in% c("Ifit3","Atp9a","Apoe","Agrn","Gpc1","Slit2"))


#### check the proteins enriched in plaque #######

plaque_meta <- turnover_meta_table %>% 
  filter(`8M and 12M plq v.s. nonplq log2FC-z` >2.5&
           `8M and 12M plq v.s. nonplq BH FDR` <0.05) %>%
  filter(`Mouse GN` %!in% c("NA")) %>% group_by(`Mouse GN`) %>% 
  slice_max(order_by = abs(`8M and 12M plq v.s. nonplq log2FC-z`), n = 1, with_ties = F)

jpeg("plaque_meta_scatter.jpeg", height = 1000, width = 1000, units = "px", res = 300)
ggplot(plaque_meta) +
  geom_point(aes(x = `WT Half-life (days)`,
                 y = `5xFAD Half-life (days)`)) +
  geom_abline(intercept = 0, slope = 1,
              color = "red", size = 1,
              alpha = 0.5) +
  #xlim(c(0,50))+
  #ylim(c(0,50)) +
  theme_bw()
dev.off()

wilcox.test(x= plaque_meta$`WT Half-life (days)`, 
            y = plaque_meta$`5xFAD Half-life (days)`)[[3]]

plaque_meta_melt <- pivot_longer(plaque_meta[,c("Mouse GN",
                                                "WT Half-life (days)",
                                                "5xFAD Half-life (days)")],
                                 cols = 2:3,
                                 names_to = "Phenotype",
                                 values_to =  "Half-life (Days)")

plaque_meta_melt <- plaque_meta_melt %>%
  filter(!is.na(`Half-life (Days)`) ) %>%
  mutate(`Log2Half-life` = log2(`Half-life (Days)`))

jpeg("plaque_meta_boxplot.jpeg", height = 1000, width = 1000, units = "px", res = 300)
ggplot(plaque_meta_melt) +
  geom_boxplot(aes(x = Phenotype, y = `Half-life (Days)`, fill = Phenotype),
               show.legend = F) +
  scale_fill_discrete()+
  theme_bw()
dev.off()




##### compare proteins in tissue specific subproteomes #########
brain_cell_type <- read_xlsx("humanBrainCellTypes_BulkRNAseq_Zhang_Neuron_2016_v1.0.xlsx", skip = 3)
human2mouse2 <- getLDS(attributes = "uniprot_gn_symbol", 
                      filters = "uniprot_gn_symbol", 
                      values = brain_cell_type$Genes, 
                      mart = human_maRt, 
                      attributesL = "uniprot_gn_symbol", 
                      martL = mouse_maRt)
names(human2mouse2) <- c("Genes","GN")
human2mouse2 <- distinct(human2mouse2, `Genes`, .keep_all = T)

brain_cell_type <- merge(brain_cell_type, human2mouse2, by = "Genes")
brain_cell_type <- merge(brain_cell_type, wt_fad_merge[,c("GN","Relative half-life change (%)")], by = "GN")
brain_cell_type[,3:8] <- log2(brain_cell_type[,3:8])


# deteremine the specificity of mRNA expression
for (i in 1:nrow(brain_cell_type)) {
  if(brain_cell_type$astrocyte[i] - max(brain_cell_type[i,4:7]) >1) {
    brain_cell_type$astrocyte_s[i] <- "Y"
  } else {
    brain_cell_type$astrocyte_s[i] <- "N"
  }
}

for (i in 1:nrow(brain_cell_type)) {
  if(brain_cell_type$neuron[i] - max(brain_cell_type[i,c(3,5:7)]) >1) {
    brain_cell_type$neuron_s[i] <- "Y"
  } else {
    brain_cell_type$neuron_s[i] <- "N"
  }
}

for (i in 1:nrow(brain_cell_type)) {
  if(brain_cell_type$oligo[i] - max(brain_cell_type[i,c(3:4,6:7)]) >1) {
    brain_cell_type$oligo_s[i] <- "Y"
  } else {
    brain_cell_type$oligo_s[i] <- "N"
  }
}

for (i in 1:nrow(brain_cell_type)) {
  if(brain_cell_type$microglia[i] - max(brain_cell_type[i,c(3:5,7)]) >1) {
    brain_cell_type$microglia_s[i] <- "Y"
  } else {
    brain_cell_type$microglia_s[i] <- "N"
  }
}

for (i in 1:nrow(brain_cell_type)) {
  if(brain_cell_type$endothelia[i] - max(brain_cell_type[i,3:6]) >1) {
    brain_cell_type$endothelia_s[i] <- "Y"
  } else {
    brain_cell_type$endothelia_s[i] <- "N"
  }
}

brain_cell_type$all_s <- "Y"
brain_cell_type_melt <- pivot_longer(brain_cell_type[,9:15], cols = 2:7,
                                     names_to = "tissue",
                                     values_to = "Y/N")
brain_cell_type_melt$tissue <- sapply(strsplit(brain_cell_type_melt$tissue,"_"),'[[',1)
brain_cell_type_melt <- filter(brain_cell_type_melt,
                               `Y/N` == "Y")
sum(brain_cell_type_melt$tissue == "all")
sum(brain_cell_type_melt$tissue == "astrocyte")
sum(brain_cell_type_melt$tissue == "neuron")
sum(brain_cell_type_melt$tissue == "oligo")
sum(brain_cell_type_melt$tissue == "microglia")
sum(brain_cell_type_melt$tissue == "endothelia")

jpeg("brain_cell_type_halflives_histogram.jpeg",
     width = 1200, height = 1000, units = "px", res = 150)
ggplot(brain_cell_type_melt) +
  geom_histogram(size = 0.5, show.legend = F,
            aes(x=  `Relative half-life change (%)`,
                 y = ..ncount..,
                 fill = tissue),alpha = 1) +
  facet_wrap(~tissue) +
  ylab("Scaled count")
dev.off()

jpeg("brain_cell_type_halflives_densityplot.jpeg",
     width = 1200, height = 1000, units = "px", res = 150)
ggplot(brain_cell_type_melt) +
  geom_density(size = 0.5, show.legend = F,
                 aes(x=  `Relative half-life change (%)`,
                     y = ..scaled..,
                     fill = tissue),alpha = 1) +
  facet_wrap(~tissue) +
  ylab("Scaled count")
dev.off()

jpeg("brain_cell_type_halflives_boxplot.jpeg",
     width = 1200, height = 1000, units = "px", res = 150)
ggplot(brain_cell_type_melt) +
  geom_boxplot(size = 0.5, show.legend = F,
                 aes(x=  tissue,
                     y = `Relative half-life change (%)`,
                     fill = tissue),alpha = 1) +
  ylab("RHC%")
dev.off()


####### compare proteins in subproteomes ##########

# protein protein target in Neuron paper
Neuron_target <- c("App","Slit2","Sfrp1","Smoc1","Htra1","Mdk", "Ntn1",
                   "Cthrc1","Ntn3","Slt3","Lsp1","C4b","Clu","Olfml3","Icam1")

Neuron_target_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% c(Neuron_target)),]

png("Neuron_target_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(Neuron_target_halflives[which(Neuron_target_halflives$WT <200),
                                           c("WT Half-life (days)",
                                             "5xFAD Half-life (days)")],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

png("Neuron_target_scatter.png", height = 1200, width = 1200,units = "px",res = 300)
ggscatterstats(Neuron_target_halflives,
               x = `WT Half-life (days)`,
               y= `5xFAD Half-life (days)`,
               label.var = GN,
               #label.expression = `log2fc of halflife (FAD vs WT)` > 0.5 
               #|`log2fc of halflife (FAD vs WT)` < -0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "white",alpha = 0),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(xlim(c(0,20)),
                      ylim(c(0,20)),
                      geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))
dev.off()





# protein up but mRNA no change in FAD v.s. WT
protein_top_up <- read.delim("../ProTransAnalysis/top_up_inconsistent_genes.txt",header = F)

protein_top_up_halflives <- wt_fad_merge[which(wt_fad_merge$GN %in% c(protein_top_up$V1)),]

jpeg("protein_top_up_boxplot.jpeg", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(protein_top_up_halflives[which(protein_top_up_halflives$WT <200),c(3,9)],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

jpeg("protein_top_up_scatter.jpeg", height = 1200, width = 1200,units = "px",res = 200)
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
jpeg("plauqe_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(melt(plaque_proteome_8m_halflives[,c(3,9)] %>% 
              filter(WT < 200,FAD < 200),
            value.name = "Half-life (days)", 
            variable.name = "mice type")) +
  geom_boxplot(aes(x = `mice type`,y= `Half-life (days)`,fill= `mice type`))
dev.off()

jpeg("plaque_scatter.png", height = 1200, width = 1200,units = "px",res = 200)
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

jpeg("overlap_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(melt(overlap_up_halflives[,c(3,9)] %>% 
              filter(WT < 200,FAD < 200),
            value.name = "Half-life (days)", 
            variable.name = "mice type")) +
  geom_boxplot(aes(x = `mice type`,y= `Half-life (days)`,fill= `mice type`))
dev.off()

jpeg("overlap_scatter.png", height = 1200, width = 1200,units = "px",res = 200)
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

jpeg("presynapse_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(presynapse_halflives[,c("WT Half-life (days)",
                                           "5xFAD Half-life (days)")],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=`Half-life (days)`)) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

jpeg("presynapse_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(presynapse_halflives,
               x = `WT Half-life (days)`,
               y= `5xFAD Half-life (days)`,
               label.var = GN,
               #label.expression = `log2fc of halflife (5xFAD vs WT)` >0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",
                                       alpha = 0, method = "lm",
                                       formula = y ~ 0 + x),
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

jpeg("M42_matrisome_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(M42_matrisome_halflives[,c("WT Half-life (days)",
                                              "5xFAD Half-life (days)")],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=log2(`Half-life (days)`))) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

jpeg("M42_matrisome_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(M42_matrisome_halflives,
               x = `WT Half-life (days)`,
               y= `5xFAD Half-life (days)`,
               label.var = GN,
               #label.expression = `log2fc of halflife (FAD vs WT)` >0.5,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",
                                       alpha = 0, method = "lm",
                                       formula = y ~ 0 + x),
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

jpeg("ABB_boxplot.png", height = 800, width = 600,units = "px",res = 200)
ggplot(data = melt(ABB_halflives[,c("WT Half-life (days)",
                                    "5xFAD Half-life (days)")],
                   variable.name = "Mice type", value.name = "Half-life (days)")
       , aes(x=`Mice type`, y=log2(`Half-life (days)`))) +
  geom_boxplot(aes(fill=`Mice type`),show.legend = F)
dev.off()

jpeg("ABB_cor.png", height = 800, width = 800,units = "px",res = 200)
ggscatterstats(ABB_halflives,
               x = `WT Half-life (days)`,
               y= `5xFAD Half-life (days)`,
               label.var = GN,
               results.subtitle = F,
               marginal = F,
               smooth.line.args = list(size = 0, color = "blue",
                                       alpha = 0, method = "lm",
                                       formula = y ~ 0 + x),
               xlab = "WT Half-lives (days)",
               ylab = "5XFAD Half-lives (days)",
               ggplot.component =
                 list(geom_abline(intercept = 0, slope = 1,
                                  color = "red", size = 1,
                                  alpha = 0.5)))                                              
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




###### human and mouse master table ######
if(FALSE) {
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
  
  human_mouse_meta_table <- merge(pro_trans_human,turnover_meta_table, by = "Mouse GN")
}
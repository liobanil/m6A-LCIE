library(png)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggprism)
library(ggrepel)
library(clusterProfiler)
library(ggbreak)
library(stringr)
library(enrichplot)

setwd("Desktop/")
#load("scatter_plot.RData")

#horizontal plot
DEGs <- read.csv("Desktop/m6A/m6A-5r/DGEA/DEG_T0vsT6.csv")
DBSs <- read.csv("/Users/gonzalal/Desktop/Poster/DBS_t6vst0.csv")
DEGs_DBSs <- read.csv("/Users/gonzalal/Desktop/LCIE/m6a-5r/DBS/t6_DEA_DBS.csv", row.names = 1)

#filter with padjs < 0.05 and |log2foldchange| > 1.5
DEGs_filtered <- DEGs %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1.5)
DBSs_filtered <- DBSs %>% filter(FDR < 0.05) %>% filter(abs(Fold) > 1)

#contrast is flipped
DBSs_filtered$Fold<-DBSs_filtered$Fold*(-1)

#horizontal bar plot
degs <- c(sum(DEGs_filtered$log2FoldChange >= 1), -1*(sum(DEGs_filtered$log2FoldChange <= -1)))
names(degs) <- c("upregulated", "downegulated")
dbss <- c(sum(DBSs_filtered$Fold  >= 1), -1*(sum(DBSs_filtered$Fold <= -1)))
names(dbss) <- c("upmethylated", "downmethylated")

df <- data.frame(NumGenes=c(degs, dbss),
                 Direction=c("L","R","L","R"),
                 Assay=c("DEG","DEG","DMS","DMS"),
                 Status=c("upregulated","downregulated","upregulated","downregulated"))

###########################################
#             Horizontal plot             #
###########################################
p1 <- ggplot(df, aes(x = Assay, y = NumGenes, fill = Direction)) +
  geom_col(position = position_stack(), width = 0.4) +
  coord_flip(clip = "off") + # Flip coordinates for horizontal bar chart
  #xlim(-70, 20) +
  scale_y_continuous(labels = abs, guide = "prism_minor") +
  labs(x = NULL, y = "number of differentially expressed/methylated genes") +
  scale_fill_manual(values = c("L" = "#4F94A6", "R" = "#B39444"),
                    labels = c("upregulated", "downregulated")) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray10") +
  theme_prism() +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "gray10", size = 20),
    axis.title = element_text(color = "gray10", size = 18), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_line(),
    plot.margin = margin(t = 5, r = 40, b = 5, l = 5),
    legend.position = "top",
    legend.direction = "horizontal", 
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.line.x = element_line(linewidth = 0.5, colour = "gray10"),
    axis.ticks.x = element_line(colour = "gray10")) + 
    geom_text(
    aes(label = abs(NumGenes)),
    hjust = ifelse(df$NumGenes < 0, -0.2, 1.2),  # Adjust text position based on count value
    vjust = 0.5,
    size = 8,
    color = "white") 
p1 
ggsave(filename = "Desktop/DEG_DBS_horizontal.png",
       plot = p1, width = 16, height = 10.33, units = "cm", dpi = 300)

###########################################
#      Tree plot for DEGs and DBs        #
###########################################
genes6 <- DEGs_DBSs$geneId[abs(DEGs_DBSs$log2FoldChange) >= 0.5] 
genes6 <- genes6[!is.na(genes6)] 
genes6 <- genes6[!duplicated(genes6)]

ego6 <- enrichGO(gene         = genes6,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP", #ALL, BP, CC, or MF
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.1,
                 readable      = TRUE)
head(ego6)

png("GO_both.png", res = 700, units = "cm", height = 12, width = 14)
dotplot(ego6, showCategory = 20, 
        font.size = 6,
        title = "DEMG - BP at 6 hours")
dev.off()

genes_sorted_6 <- both$log2FoldChange[!duplicated(genes6)]
names(genes_sorted_6) <- both$geneId[!duplicated(genes6)]
genes_sorted_6 <- na.omit(genes_sorted_6)
genes_sorted_6 <- sort(genes_sorted_6, decreasing = TRUE)


gse_6 <- gseGO(geneList=genes_sorted_6, 
               ont ="BP", 
               keyType = "ENTREZID",
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "none")

head(gse_6)

neural_terms <- gse_6@result %>%
  filter(grepl("neuro|synapse|axon|brain|neural", Description, ignore.case=TRUE))


gse_6@result <- neural_terms

png(paste0("gse_6_",enrich,"_FC_1.png"), res = 700, units = "cm", height = 10, width = 18)
print({
  dotplot(gse_6, showCategory=, split=".sign", font.size = 10) + theme(legend.position = "right") + facet_grid(.~.sign)  + scale_y_discrete(labels=function(egoBP) str_wrap(egoBP, width=60))
})
dev.off()

ego.d <- pairwise_termsim(ego6)

png(paste("Desktop/DMEG_ego_6_treeplot_BP.png", sep = ""), res = 700, units = "cm", height = 18, width = 35)
treeplot(ego.d, cluster.params = list(label_words_n = 2, method = "average"), hilight.params = list(hilight=TRUE, align="both")) +
  theme(legend.position="right")
dev.off()

###########################################
#     Scatter plot: DBSs and DGEs         #
###########################################

#####
DBS_annot <- read.csv("/Users/gonzalal/Desktop/Poster/m6a_DBA_t6vst0_annot.csv", row.names = 1)
DBS_annot$Fold<-DBS_annot$Fold*(-1) #contrast flipped
colnames(DEGs)[1] <- "SYMBOL"

#idx_deg <- DEGs$SYMBOL %in% DBS_annot$SYMBOL
#idx_dbs <- DBS_annot$SYMBOL %in% DEGs$SYMBOL

DEGs_DBSs <- read.csv("/Users/gonzalal/Desktop/LCIE/m6a-5r/DBS/t6_DEA_DBS.csv", row.names = 1)
DEGs_DBSs <- DEGs_DBSs[,c(1,4,6)]

# Set limits for x and y axes
x_limit <- c(-2.2, 4.5)
y_limit <- c(-2.5, 5)

#genes of interest #Mct4 (Slc16a3) is overexpressed
#genes_interest <- c("Bdnf", "Ythcd1","Ngf", "Kctd8", "Kctd4", "Vegfa")
neurotrofinas <- c("Bdnf","Ngf",  "Vegfa")
k_channles <- DEGs_DBSs$SYMBOL[grep("Kc", DEGs_DBSs$SYMBOL)]
ca_channles <- DEGs_DBSs$SYMBOL[grep("Cacna", DEGs_DBSs$SYMBOL)]
m6a_enz <- "Ythdc1"
plasticity <- c("Fosb", "Adcy8")
genes_interest <- c(neurotrofinas,k_channles,ca_channles,m6a_enz,plasticity)

#add color column to highlight relevant genes
DEGs_DBSs$color <- ""

DEGs_DBSs$color <- "none"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% neurotrofinas] <- "B"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% k_channles] <- "G"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% ca_channles] <- "P"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% m6a_enz] <- "C"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% plasticity] <- "Y"
DEGs_DBSs$color <- factor(DEGs_DBSs$color, levels = c("Y", "C", "P", "B", "G", "none"))

#subset ego.d
#####
top30go <- as.data.frame(ego.d)[1:30,]
top30go$color <- "#D3C781"
top30go$color[grep("guidance", top30go$Description)] <- "#E1A4CA"
top30go$color[grep("eye|ERK1|morphogenesis|visual|sensory", top30go$Description)] <- "#8CD3D7"
top30go$color[grep("protein modification", top30go$Description)] <- "#ACC0E3"
top30go$color[grep("muscle cell", top30go$Description)] <- "#9DCD92"


yellow_genes <- unique(unlist(stringr::str_split(paste(top30go$geneID[top30go$color=="#D3C781"], 
                                                collapse = "/"),pattern = "/")))
green_genes <- unique(unlist(stringr::str_split(paste(top30go$geneID[top30go$color=="#9DCD92"], 
                                               collapse = "/"),pattern = "/")))
pink_genes <- unique(unlist(stringr::str_split(paste(top30go$geneID[top30go$color=="#E1A4CA"], 
                                              collapse = "/"),pattern = "/")))
cyan_genes <- unique(unlist(stringr::str_split(paste(top30go$geneID[top30go$color=="#8CD3D7"], 
                                              collapse = "/"),pattern = "/")))
blue_genes <- unique(unlist(stringr::str_split(paste(top30go$geneID[top30go$color=="#ACC0E3"], 
                                              collapse = "/"),pattern = "/")))

DEGs_DBSs$color <- "none"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% blue_genes] <- "B"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% green_genes] <- "G"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% pink_genes] <- "P"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% cyan_genes] <- "C"
DEGs_DBSs$color[DEGs_DBSs$SYMBOL %in% yellow_genes] <- "Y"
DEGs_DBSs$color <- factor(DEGs_DBSs$color, levels = c("Y", "C", "P", "B", "G", "none"))

#####
cidx <- DEGs_DBSs$SYMBOL %in% genes_interest
#DEGs_DBSs$color[cidx] <- "Of interest"


#set shape and reset coordinates if out of range
# out of range on x
#sum(DEGs_DBSs$log2FoldChange > 6 | DEGs_DBSs$log2FoldChange < -3) - none out of range
# out of range in y
# sum(DEGs_DBSs$Fold > 5 | DEGs_DBSs$Fold < -2.5) - 6 out of range
ridxp <- DEGs_DBSs$Fold > 5
DEGs_DBSs$Fold[ridxp] <- 5
ridxn <- DEGs_DBSs$Fold < -2.5
DEGs_DBSs$Fold[ridxn] <- -2.5

#add column for shape
DEGs_DBSs$shape <- 1
DEGs_DBSs$shape[DEGs_DBSs$color == "none"] <- 4
DEGs_DBSs$shape[ridxp] <- 2
DEGs_DBSs$shape[ridxn] <- 3
DEGs_DBSs$shape <- as.factor(DEGs_DBSs$shape)


#scatter plot
p <- ggplot(DEGs_DBSs, aes(x = log2FoldChange, y = Fold, color = color)) +
  geom_point(size=3,aes(shape=shape)) + # Use shape 1 (hollow circle) for points
  xlim(x_limit) + ylim(y_limit) +    # Set x and y axis limits
  theme_minimal() +  # Optional: Use a minimal theme
  labs(x = "Expression log2FoldChange", y = "Methylation FoldChange") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16)) +
  labs(x = "Expression (log2 fold change)",
       y = "Methylation (log2 fold change)") +
  scale_color_manual(values = c("G" = "#9DCD92", "C" = "#FF746C", 
                    "Y" = "#D3C781", "P" = "#E1A4CA", "B" = "#ACC0E3","none" = "black"), 
                    labels = c("G" = "Potassium channels", "C" = "Ythdc1", 
                               "Y" = "Plasticity related genes", "P" = "Calcium channels", 
                               "B" = "Neurotrophins", "none" = "Other genes")) +
  scale_shape_manual("Shape code",values = c("1"=16, "2"=24, "3"=25, "4"=1)) +
  theme(axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 16),  # Increase legend title size
        legend.text = element_text(size = 14),  # Increase legend text size
        legend.key.size = unit(1.5, "lines"),   # Increase the size of legend keys
        legend.position = "bottom") +
  guides(shape = FALSE) +
  geom_label_repel(data = subset(DEGs_DBSs, SYMBOL %in% genes_interest), 
                   aes(label = SYMBOL), 
                   size = 3.5, 
                   color = "#FF6961",
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines"),
                   segment.color = "#FF6961",
                   segment.size = 0.5,  
                   nudge_y = -DEGs_DBSs[cidx,]$log2FoldChange*0.5,
                   nudge_x = 0.5*DEGs_DBSs[cidx,]$Fold, 
                   force_pull = 1, 
                   xlim = c(-3,6), 
                   ylim = c(-2.5,5))

  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", linewidth = 1.2)

   

# Display the plot
p 
ggsave(filename = "Desktop/DEG_DBS_scatter_plot_labels.png",
       plot = p, width = 25, height = 15, units = "cm", dpi = 300)

#save
save.image(file = "scatter_plot.RData")
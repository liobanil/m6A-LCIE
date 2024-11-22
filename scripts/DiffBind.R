library(GenomicRanges)
library(DiffBind)
library(Guitar)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(GenomicFeatures)
library(rtracklayer)
library(EnsDb.Mmusculus.v79)

#####################################################################
#                              ChIPQC                               #
#####################################################################
require(ChIPQC)

#ceate chipqc object
chipObj <- ChIPQC("sampleinfo.csv", annotation="mm10") 

#make report
ChIPQCreport(chipObj, reportName="ChIP QC Report", reportFolder="DiffBind/ChIPQCreport")

#higher SSD is more indicative of better enrichment
#FRIP: percentage of reads that overlap called peaks. Indication of how "enriched" the sample is, and can be considered a "signal-tonoise" measure of what proportion of the library consists of fragments from binding sites vs. background reads.
#RiP% values for ChIPs around 5% or higher generally reflect successful enrichment.
#####################################################################
#                             DiffBind                              #
#####################################################################

require(DiffBind)

#load samples
print("Reading samplesheet into DBA object...")
samples <- read.csv("/home/gonzalal/scratch/m6a-5r/DiffBind/sampleinfo.csv", header = TRUE, row.names = 1)
samples$SampleID <- make.names(samples$SampleID)
samples$ControlID <- make.names(samples$ControlID)

#create dba object
m6a <- dba(sampleSheet = samples)

print("Generating Heatmaps...")
png("Correlation_hetmap.png", )
plot(m6a)
dev.off()

#dba coutns and keep peaks only present in at least 3 replicates
print("Computing counts...")
m6a <- dba.count(m6a,  bUseSummarizeOverlaps = FALSE, minOverlap = 2, score="DBA_SCORE_READS")

#save count matrix
print("Creating count files...")
counts_m6a <- dba.peakset(m6a, bRetrieve=TRUE, writeFile="DiffBind/Counts_m6a.txt")

#normalize
m6a <- dba.normalize(m6a, method = DBA_DESEQ2, normalize = DBA_NORM_RLE)

#Model design and contrast
print("Making contrast and adding blacklists...")
#m6a <- dba.contrast(m6a, reorderMeta = list(Treatment="0"))
m6a <- dba.contrast(m6a, contrast=c("Treatment", "6", "0"))
m6a <- dba.contrast(m6a, contrast=c("Treatment", "24", "0"))

#Differential Binding Analysis
print("Binding analysis...")
m6a <- dba.analyze(m6a, method = DBA_DESEQ2, bBlacklist = DBA_BLACKLIST_MM10, bGreylist = FALSE)

#Retrieve the differentially methylated sites
m6a_t0_t6 <- dba.report(m6a, contrast = 1) 
m6a_t0_t24 <- dba.report(m6a, contrast = 2) 

#write csv with DMSs
m6a_t0_t6_df <- as.data.frame(m6a_t0_t6)
m6a_t0_t24_df <- as.data.frame(m6a_t0_t24)

write.csv(m6a_t6_DBS_annot_df,"/home/gonzalal/scratch/m6a-5r/DiffBind/DMSs/DMS_t6_vs_t0.csv")
write.csv(m6a_t24_DBS_annot_df,"/home/gonzalal/scratch/m6a-5r/DiffBind/DMSs/DMS_t24_vs_t0.csv")

#heatmap of DMSs
print("Making heatmaps of DMSs...")
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

png("m6A_DBA_Heatmap_t6vst0.png", res = 700 ,width=15.23, height=20.78, units = "cm")
dba.plotHeatmap(m6a, contrast=1, correlations=FALSE, scale="row", colScheme=hmap)
dev.off()

png("m6A_DBA_Heatmap_t24vst0.png", res = 700 ,width=15.23, height=20.78, units = "cm")
dba.plotHeatmap(m6a, contrast=2, correlations=FALSE, scale="row", colScheme=hmap)
dev.off()


#####################################################################
#                        Peak annotation                            #
#####################################################################
print("Annotating differentially binding sites...")

#load annotation..
print("Loading annotation...")
annoData <- toGRanges(EnsDb.Mmusculus.v79, feature = "gene")
chr <- read.table("/home/gonzalal/scratch/m6a-5r/mm39/chr_equivalence.txt")

#change chr names to match annotation
chr_m6a <- chr$V2[match(as.character(seqnames(m6a)), chr$V1)]
chr_t0_t6 <- chr$V2[match(as.character(seqnames(m6a_t0_t6)), chr$V1)]
chr_t0_t24 <- chr$V2[match(as.character(seqnames(m6a_t0_t24)), chr$V1)]
all <-  GRanges(chr_m6a, ranges(m6a))
t0_t6 <- GRanges(chr_t0_t6, ranges(m6a_t0_t6))
t0_t24 <- GRanges(chr_t0_t24, ranges(m6a_t0_t24))

#annotate overlaps
print("Annotate overlaping peak regions within genes...")
#all
overlaps.all <- annotatePeakInBatch(all, AnnotationData=annoData, output="overlapping", maxgap=5000L)
overlaps.all$gene_name <- annoData$gene_name[match(overlaps.all$feature, names(annoData))]
#t0_t6
overlaps.t0_t6 <- annotatePeakInBatch(t0_t6, AnnotationData=annoData, output="overlapping", maxgap=5000L)
overlaps.t0_t6$gene_name <- annoData$gene_name[match(overlaps.t0_t6$feature, names(annoData))]
#t0_t24
overlaps.t0_t24 <- annotatePeakInBatch(t0_t24, AnnotationData=annoData, output="overlapping", maxgap=5000L)
overlaps.t0_t24$gene_name <- annoData$gene_name[match(overlaps.t0_t24$feature, names(annoData))]

#write anoteted peaks files
print("Creating annotated peak files...")
rtracklayer::export.bed(overlaps.all, "DiffBind/PeakAnnotation_all.bed")
rtracklayer::export.bed(overlaps.t0_t6, "/home/gonzalal/scratch/m6a-5r/DiffBind/DMSs/PeakAnnotation_t0_t6.bed")
rtracklayer::export.bed(overlaps.t0_t24, "/home/gonzalal/scratch/m6a-5r/DiffBind/DMSs/PeakAnnotation_t0_t24.bed")



#with ChIPseeker
m6a_t6_DBS_df <- as.data.frame(m6a_t0_t6)
#change seqlevels
chr <- read.table("/home/gonzalal/scratch/m6a-5r/DiffBind/chr_eq.txt")
m6a_t6_DBS_df$seqnames <- chr$V2[match(m6a_t6_DBS_df$seqnames, chr$V1)]

#select abs(foldchange) > 2 and FDR < 0.01
m6a_t6_DBS_selected_df <-  m6a_t6_DBS_df[m6a_t6_DBS_df$FDR < 0.01 & abs(m6a_t6_DBS_df$Fold) > 2,]
m6a_t6_DBS_selected <- makeGRangesFromDataFrame(m6a_t6_DBS_selected_df, keep.extra.columns=TRUE)

#annotate the sites
m6a_t6_DBS_selected_annot <- annotatePeak(m6a_t6_DBS_selected, tssRegion=c(-3000,3000), TxDb=TxDb.Mmusculus.UCSC.mm39.refGene, annoDb="org.Mm.eg.db")

#make df
m6a_t6_DBS_selected_annot_df <- as.data.frame(m6a_t6_DBS_selected_annot)
write.csv(m6a_t6_DBS_annot_df,"DMS_t6_vs_t0_annot.csv")

#makeGRanges again
m6a_t6_DBS_UCSC <- makeGRangesFromDataFrame(m6a_t6_DBS_df)

#copy metadata
m6a_t6_DBS_UCSC$Conc <- m6a_t6_DBS$Con
m6a_t6_DBS_UCSC$Conc_0 <- m6a_t6_DBS$Conc_0
m6a_t6_DBS_UCSC$Conc_6 <- m6a_t6_DBS$Conc_6
m6a_t6_DBS_UCSC$Fold <- m6a_t6_DBS$Fold
m6a_t6_DBS_UCSC$"p-value" <- m6a_t6_DBS$"p-value"
m6a_t6_DBS_UCSC$FDR <- m6a_t6_DBS$FDR

#annotate peaks with ChIPseeker
m6a_t6_DBS_annot <- annotatePeak(m6a_t6_DBS_UCSC, tssRegion=c(-3000,3000), TxDb=TxDb.Mmusculus.UCSC.mm39.refGene, annoDb="org.Mm.eg.db")

#make df
m6a_t6_DBS_annot_df <- as.data.frame(m6a_t6_DBS_annot)
write.csv(m6a_t6_DBS_annot_df,"m6a_t6_DBS_annot_df.csv")


#t24
m6a <- dba.contrast(m6a, contrast=c("Treatment", "24", "0"), reorderMeta=list(Treatment="0"))

m6a <- dba.analyze(m6a, method=DBA_ALL_METHODS, bBlacklist=DBA_BLACKLIST_MM10)

m6a_t24_DBS <- dba.report(m6a, contrast=2)

m6a_t24_DBS_df <- as.data.frame(m6a_t24_DBS)

#change seqlevels
m6a_t24_DBS_df$seqnames <- chr$V2[match(m6a_t24_DBS_df$seqnames, chr$V1)]

m6a_t24_DBS_selected_df <-  m6a_t24_DBS_df[m6a_t24_DBS_df$FDR < 0.01 & abs(m6a_t24_DBS_df$Fold) > 2,]
m6a_t24_DBS_selected <- makeGRangesFromDataFrame(m6a_t24_DBS_selected_df, keep.extra.columns=TRUE)

m6a_t24_DBS_selected_annot <- annotatePeak(m6a_t24_DBS_selected, tssRegion=c(-3000,3000), TxDb=TxDb.Mmusculus.UCSC.mm39.refGene, annoDb="org.Mm.eg.db")

m6a_t24_DBS_selected_annot_df <- as.data.frame(m6a_t24_DBS_selected_annot)
write.csv(m6a_t24_DBS_selected_annot_df,"m6a_t24_DBS_annot_df.csv")

#####################################################################
#                           Guitar Plots                            #
#####################################################################
print("Making gene feature plots (guitar plots)...")

#make .narrowPeak into .bed
stBedFiles <- list(file.path("Guitar/Guitar_tmp.bed"))
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_narrowPeak <- list(import("../MACS/Ctl_1/Ctl_1_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak))

ctl <- c("../MACS/Ctl_1/Ctl_1_peaks.narrowPeak", "../MACS/Ctl_2/Ctl_2_peaks.narrowPeak", "../MACS/Ctl_3/Ctl_3_peaks.narrowPeak", "../MACS/Ctl_4/Ctl_4_peaks.narrowPeak", "../MACS/Ctl_5/Ctl_5_peaks.narrowPeak")

t6 <- c("../MACS/t6_1/t6_1_peaks.narrowPeak", "../MACS/t6_2/t6_2_peaks.narrowPeak", "../MACS/t6_3/t6_3_peaks.narrowPeak","../MACS/t6_4/t6_4_peaks.narrowPeak", "../MACS/t6_5/t6_5_peaks.narrowPeak")

t24 <- c("../MACS/t24_1/t24_1_peaks.narrowPeak", "../MACS/t24_2/t24_2_peaks.narrowPeak", "../MACS/t24_3/t24_3_peaks.narrowPeak", "../MACS/t24_4/t24_4_peaks.narrowPeak", "../MACS/t24_5/t24_5_peaks.narrowPeak")

ctl <- lapply(ctl, import, extraCols=extraCols_narrowPeak, format="BED")
ctl <- GRangesList(ctl)
ctl_red <- reduce(ctl, min.gapwidth=31)
ctl <-  do.call("c", ctl_red)
ctl <- reduce(ctl, min.gapwidth=31)
chr_names <- chr$V2[match(as.character(seqnames(ctl)),chr$V1)]
ctl.2 <- GRanges(chr_names, ranges(ctl), strand(ctl))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
ctl.2 <- keepStandardChromosomes(ctl.2, pruning.mode="coarse")
ctl.2.l <- GRangesList(ctl.2)
GuitarPlot(txTxdb=txdb, stGRangeLists=ctl.2.l, miscOutFilePrefix = "Guitar Plot Ctl_1")

t6 <- lapply(t6, import, extraCols=extraCols_narrowPeak, format="BED")
t6 <- GRangesList(t6)
t6_red <- reduce(t6, min.gapwidth=31)
t6 <-  do.call("c", t6_red)
t6 <- reduce(t6, min.gapwidth=31)
chr_names <- chr$V2[match(as.character(seqnames(t6)),chr$V1)]
t6.2 <- GRanges(chr_names, ranges(t6), strand(t6))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
t6.2 <- keepStandardChromosomes(t6.2, pruning.mode="coarse")
t6.2.l <- GRangesList(t6.2)
GuitarPlot(txTxdb=txdb, stGRangeLists=t6.2.l, miscOutFilePrefix = "Guitar Plot t6")

t24 <- lapply(t24, import, extraCols=extraCols_narrowPeak, format="BED")
t24 <- GRangesList(t24)
t24_red <- reduce(t24, min.gapwidth=31)
t24 <-  do.call("c", t24_red)
t24 <- reduce(t24, min.gapwidth=31)
chr_names <- chr$V2[match(as.character(seqnames(t24)),chr$V1)]
t24.2 <- GRanges(chr_names, ranges(t24), strand(t24))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
t24.2 <- keepStandardChromosomes(t24.2, pruning.mode="coarse")
t24.2.l <- GRangesList(t24.2)
GuitarPlot(txTxdb=txdb, stGRangeLists=t24.2.l, miscOutFilePrefix = "Guitar Plot t24")









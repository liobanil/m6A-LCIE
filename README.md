# Lactate influences the m6A methylation and expression of genes involved in neuronal function
Author: *AdLLG*

Introduction
------
Lactate is both a fuel and a signaling molecule that plays important roles in neuronal function and development. In the past decade, the post-transcriptional mark commonly referred to as m6A, a methylation in the position N6 of the adenosine of RNA, has been recognized to be involved in important processes such as development and cognitive functions, such as learning and memory. It has not been completely described how lactate executes its function. Given the overlap between the m6A mark and lactate roles, an interplay between them can help elucidate better how lactate sustains neuronal function.

Experimental design
-----
Neuronal primary cultures were treated with lactate (10 mM) for 6 and 24 hours. Then RNA was extracted and ployA tail was used to select mRNA. Subsequently, immunoprecipitation of m6A-marked RNA molecules was performed and libraries were made. The samples were sequenced using Illumina 150 bp paired-end reads for both libraries. 

Protocols: 
MeRIP: [EpiQuik™ CUT&RUN m6A RNA Enrichment (MeRIP) Kit](https://www.epigentek.com/docs/P-9018.pdf)

Library preparation: [Takara SMARTer Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian](https://www.takarabio.com/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA/SMARTer%20Stranded%20Total%20RNA-Seq%20Kit%20v2%20-%20Pico%20Input%20Mammalian%20User%20Manual_050619.pdf)

Experimental design summary:
5 replicates per condition (control, 6hrs, and 24 hrs).

**Bioinformatics analysis:**



Directories and content 
-----

1. original_data

    Original sequencing data along with FastQC results. 

2. trimmed_data

    Illumina Universal adapter removal using cutadapt along with FastQC results.

3. STAR

    .bam obtained with STAR

4. HISAT2

    .bam obtained with Hisat2

5. MACS2

    Results from running macs2

6. HOMER

    Results from de novo motif discovery

7. DiffBind

    QC results performed with ChIPQC and plots and results from R script using DiffBinf
    
8. scripts

    Bash and R scripts used in this pipeline
    
9. stdout

    Relevant stdout and err messages
    
10. mm39

    Genome and annotation used to run the analysis


Pipeline
-----

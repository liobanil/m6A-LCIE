# m6A-LCIE
Author: *AdLLG*

Description
------
Scripts and files to run the pipeline to analyze differential binding on meRIP seq data. 

Experimental design
-----
Neuronal primary cultures were treated with lactate (10 mM) for 6 and 24 hours. Then RNA was extracted and ployA tail was used to select the mRNA. The m6A immunoprecipitation was performed and libraries were made. The samples were sequenced using Illumina 150 bp paired-end reads for both libraries (Core Lab or NovoGene). 

Protocols: 
MeRIP: [EpiQuik™ CUT&RUN m6A RNA Enrichment (MeRIP) Kit](https://www.epigentek.com/docs/P-9018.pdf)

Library preparation: [Takara SMARTer Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian](https://www.takarabio.com/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA/SMARTer%20Stranded%20Total%20RNA-Seq%20Kit%20v2%20-%20Pico%20Input%20Mammalian%20User%20Manual_050619.pdf)

Experimental design for the first library:
3 replicates per condition (control, 6hrs, and 24 hrs).
<img width="905" alt="exp_desing" style="width:400px;" src="https://user-images.githubusercontent.com/54646526/194641784-d18dda3b-521d-41ab-babd-9d014677954b.png">

Experimental design for the second library:
5 replicates per condition (control, 6hrs, and 24 hrs).

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
# Lactate influences the m6A methylation and expression of genes involved in neuronal function
Author: *AdLLG* :shipit:

## Introduction
Lactate is both a fuel and a signaling molecule that plays important roles in neuronal function and development. In the past decade, the post-transcriptional mark commonly referred to as m6A, a methylation in the position N6 of the adenosine of RNA, has been recognized to be involved in important processes such as development and cognitive processes, such as learning and memory. In this project, we aim to explore the relationship between lactate and the m6A mark.

The purpose of this GitHub is to document the pipeline used to analyze and generate the results for the m6A data as well as the contents of the directory.

Experimental design
-----

![Experimental design and methodology description.](/img/Methods_Fig.png)


Cortical neuronal primary cultures were treated with lactate (10 mM) for 6 and 24 hours. Total RNA was extracted followed by polyA selection to enrich for mRNA. Subsequently, immunoprecipitation of m6A-marked mRNA molecules was performed and sequencing libraries were made. The samples were sequenced using Illumina 150 bp paired-end reads for both libraries. 

Experimental design summary:
control, 6 hours lactate-treated, and 24 hours lactate-treated, each group counts with 5 replicates. 


Protocols: 
MeRIP: [EpiQuik™ CUT&RUN m6A RNA Enrichment (MeRIP) Kit](https://www.epigentek.com/docs/P-9018.pdf)

Library preparation: [Takara SMARTer Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian](https://www.takarabio.com/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA/SMARTer%20Stranded%20Total%20RNA-Seq%20Kit%20v2%20-%20Pico%20Input%20Mammalian%20User%20Manual_050619.pdf)

## Bioinformatic Analysis

The data analysis was performed on Ibex. The QC of the raw data, adapter removal (trimming), mapping of processed data, and peak calling were performed by running snakemake. The folder _rules_ contains `.smk` for each step. For running, the rules were called with `.sh` which invokes snakemake in the following way: `snakemake --snakefile snakefile --use-conda`. The _snakefile_ is available inside the _scripts_ folder. 

Differential gene expression analysis (results available in folder _DGEA_) and differential methylation analysis (results available in folder _DiffBind_) were performed with R scripts, available in _scripts_ folder. 

The _conda_ folder contains the `.yaml` files for the mapping and trimming steps. The logs for each step were stored in the _logs_ folder (empty since its not relevant to store those files here). _mm39_ and _mm10_ contained the annotation used for some analysis steps, both versions are used since some R packages still depend on mm10. _config_ folder has the `.yaml` file for snakemake.

The _qPCR_ and _MEA_ folders contain the R scripts to analyze the qPCR and neural activity recordings, respectively. These folders are NOT on Ibex.

Directories and content 
-----
> [!NOTE]
> Raw, trimmed, and MACS data is unavailable here due to size restrictions.

The following tree describes what each folder contains, recreating the working directory on Ibex. 

```
m6A-LCIE 
  ├── raw_data
  ├── trimmed_data
  ├── mapped_data
  ├── MACS
  ├── conda
  ├── config
  ├── scripts
  ├── logs
  ├── mm39
  ├── mm10
  ├── rules
  ├── scripts
  ├── DiffBind
  ├── DGEA
  ├── qPCR
  └── MEA
```

**How to run the pipeline:**

```

```

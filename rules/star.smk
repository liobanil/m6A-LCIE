rule star_index:
    input:
        fasta=config["mm39"]["fna"]
    output:
        directory("mm39/star_index_50")
    threads: 1
    resources:
        mem_mb=168432
    params:
        extra="--sjdbGTFfile mm39/annotation/ncbi-genomes-2023-06-01/GCF_000001635.27_GRCm39_genomic.gtf --sjdbOverhang 50 --limitGenomeGenerateRAM 168422653312",
    log:
        "logs/star_index/star_index.log"
    wrapper:
        "v2.2.1/bio/star/index"

rule loadandkeep_index:
    input:
        index=config["mm39"]["index"]
    output:
        done=touch("loading.done")
    shell:
        """
        STAR --genomeLoad LoadAndExit --genomeDir {input.index}
        touch {output.done}
        """

rule star_mapping:
    input:
        fq1 = config["trimmed_data"] + "{sample}/{sample}_1.fq.gz",
        fq2 = config["trimmed_data"] + "{sample}/{sample}_2.fq.gz",
        idx = config["mm39"]["index"]
    params:
        extra=" --outSAMunmapped Within \
                --outSAMattributes Standard \
                --quantMode GeneCounts \
                --outSAMtype BAM SortedByCoordinate"
    output:
        aln = config["mapped_data"] + "{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        reads_per_gene = config["mapped_data"] + "{sample}/{sample}_ReadsPerGene.out.tab",
        sj = config["mapped_data"] + "{sample}/{sample}_SJ.out.tab",
        log = config["mapped_data"] + "{sample}/{sample}_Log.out",
        log_progress = config["mapped_data"] + "{sample}/{sample}_Log.progress.out",
        log_final = config["mapped_data"] + "{sample}/{sample}_Log.final.out"
    log:
        "logs/star_mapping/star_mapping_{sample}.log"
    threads: 20
    wrapper:
        "v2.2.1/bio/star/align"

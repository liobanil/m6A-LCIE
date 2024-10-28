rule trimm_raw_data:
    input:
        fq1=config["raw_data"] + "{sample}/{sample}_DKDL230006174-1A_HG3H7DSX7_L3_1.fq.gz",
        fq2=config["raw_data"] + "{sample}/{sample}_DKDL230006174-1A_HG3H7DSX7_L3_2.fq.gz"
    output:
        fastq1=config["trimmed_data"] +"{sample}/{sample}_1.fq.gz",
        fastq2=config["trimmed_data"] +"{sample}/{sample}_2.fq.gz",
        qc=config["trimmed_data"] +"{sample}/{sample}.qc.txt"
    conda: "conda/trimming-env.yaml"
    params:
        adap_a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adap_A= "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    log:
        "logs/trimm_raw_data/{sample}.log"
    threads: 10
    shell:
        "cutadapt -a {params.adap.a} -A {params.adap_A} -o {output.fastq1} -p {output.fastq2} {input.fq1} {input.fq2} --minimum-length 50 -q 20"

rule fastqc_trimm_data:
    input:
        config["trimmed_data"] + "{sample}/{sample}_{read_no}.fq.gz"
    output:
        html=config["trimmed_data"] + "qc/fastqc/{sample}_{read_no}.html",
        zip=config["trimmed_data"] + "qc/fastqc/{sample}_{read_no}_fastqc.zip"
    params:
        extra="--quiet"
    log:
        "logs/fastqc_trimm_data/{sample}_{read_no}.log"
    threads: 10
    resources:
        mem_mb=1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule multiqc_trimm_data:
    input:
        expand(config["trimmed_data"] + "qc/fastqc/{sample}_{read_no}.html", sample=SAMPLES, read_no=['1','2'])
    output:
        config["trimmed_data"] + "qc/multiqc.html"
    log:
        "logs/multiqc_trimm_data/multiqc.log"
    wrapper:
        "v2.2.1/bio/multiqc"

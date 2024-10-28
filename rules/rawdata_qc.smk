rule fastqc_raw_data:
    input:
        config["raw_data"] + "{sample}/{sample}_DKDL230006174-1A_HG3H7DSX7_L3_{read_no}.fq.gz"
    output:
        html=config["raw_data"] + "qc/fastqc/{sample}_{read_no}.html",
        zip=config["raw_data"] + "qc/fastqc/{sample}_{read_no}_fastqc.zip"
    params:
        extra="--quiet"
    log:
        "logs/fastqc_raw_data/{sample}_{read_no}.log"
    threads: 10
    resources:
        mem_mb=1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule multiqc_raw_data:
    input:
        expand(config["raw_data"] + "qc/fastqc/{sample}_{read_no}.html", sample=SAMPLES, read_no=['1','2'])
    output:
        config["raw_data"] + "qc/multiqc.html"
    log:
        "logs/multiqc_raw_data/multiqc.log"
    wrapper:
        "v2.2.1/bio/multiqc"

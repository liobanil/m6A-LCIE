rule callpeak:
    input:
        treatment=config["mapped_data"] + "IP_{sample}/IP_{sample}_Aligned.sortedByCoord.out.bam",   # ip samples
        control=config["mapped_data"] + "Input_{sample}/Input_{sample}_Aligned.sortedByCoord.out.bam"      # input samples
    output:
        multiext("MACS/{sample}/{sample}",
                 "_peaks.xls",
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 "_summits.bed",
                 "_peaks.narrowPeak"
                 )
    log:
        "logs/callpeak/{sample}_callpeak.log"
    params:
        "-f BAM -g 2.5e+9 --nomodel --bdg"
    wrapper:
        "v2.2.1/bio/macs2/callpeak"

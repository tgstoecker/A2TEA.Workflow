ruleorder: Trimmomatic_PE > Trimmomatic_SE


def get_trimmed_PE_output(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"trimmed/{s.fq1}", f"trimmed/unpaired/unpaired_{s.fq1}", f"trimmed/{s.fq2}", f"trimmed/unpaired/unpaired_{s.fq2}"  ]


def get_trimmed_SE_output(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"trimmed/{s.fq1}"  ]



#this will create the directories needed for the trimmomatic rules
TRIM_DIR = "trimmed/unpaired/"
if not os.path.exists(TRIM_DIR):
    os.makedirs(TRIM_DIR)

rule Trimmomatic_PE:
    input:
        get_PE_fastqs,
    output:
        touch("checks/trimmed/{sample}_{unit}_PE.check"),
    log:
        "logs/trimmomatic/PE/{sample}_{unit}.log"
    params:
        illuminaclip=config["trim_illuminaclip"],
        trimmer=config["trim_options"],
        compression_level="-9",
        out = get_trimmed_PE_output,
    threads:
        config["threads_trimmomatic"]
    shell:
        "trimmomatic PE -threads {threads} {input} {params.out} {params.illuminaclip} {params.trimmer} 2> {log}"


rule Trimmomatic_SE:
    input:
        get_SE_fastqs,
    output:
        touch("checks/trimmed/{sample}_{unit}_SE.check"),
    log:
        "logs/trimmomatic/SE/{sample}_{unit}.log"
    params:
        illuminaclip=config["trim_illuminaclip"],
        trimmer=config["trim_options"],
        compression_level="-9",
        out = get_trimmed_SE_output,
    threads:
        config["threads_trimmomatic"]
    shell:
        "trimmomatic SE -threads {threads} {input} {params.out} {params.illuminaclip} {params.trimmer} 2> {log}"


rule Trimmomatic_cleanup:
    input:
        expand("checks/trimmed/{sample.sample}_{sample.unit}_SE.check", sample=SE_samples.itertuples()),
        expand("checks/trimmed/{sample.sample}_{sample.unit}_PE.check", sample=PE_samples.itertuples()),
    output:
        touch("checks/trimmed/trim_cleanup.check"),
    shell:
        "rm -r trimmed/unpaired/"

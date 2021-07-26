#this will create the directories needed for the trimmomatic rules
TRIM_DIR = "trimmed/unpaired/"
if not os.path.exists(TRIM_DIR):
    os.makedirs(TRIM_DIR)


def get_trimmed_PE_output(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"trimmed/{s.fq1}", f"trimmed/unpaired/unpaired_{s.fq1}", f"trimmed/{s.fq2}", f"trimmed/unpaired/unpaired_{s.fq2}"  ]


def get_trimmed_SE_output(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"trimmed/{s.fq1}"  ]


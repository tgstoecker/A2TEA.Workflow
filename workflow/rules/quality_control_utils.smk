#this will create (if needed) the directories for the trimmomatic rules
FQC_RAW_DIR = "fastqc/raw/"
if not os.path.exists(FQC_RAW_DIR):
    os.makedirs(FQC_RAW_DIR)

FQC_TRIM_DIR = "fastqc/trimmed/"
if not os.path.exists(FQC_TRIM_DIR):
    os.makedirs(FQC_TRIM_DIR)

# if there is at least 1 species for which a genomic fasta was supplied -> classic alignment with STAR
if len(GEN_FASTA_SPECIES) != 0:


    def get_paired_info(wildcards):
        """Get raw FASTQ files from unit sheet."""
        if not is_single_end(wildcards.sample, wildcards.unit):
            return [ f"-p" ]
        else:
            return [ f"" ]

    def getCountsForSpecies(wildcards):
        counts = list()
        # only the actual counts (ending with .txt) are considered in the diectory
        # also the output is sorted - when sticking to control and treatment everything is as it should
        # control 1-4, treatment 1-4
        for c in sorted(os.listdir("featureCounts/"+wildcards.species+"/gene_level/")):
            if c.endswith(".txt"):
                counts.append(os.path.join("featureCounts/",wildcards.species,"gene_level/",c))
        return counts


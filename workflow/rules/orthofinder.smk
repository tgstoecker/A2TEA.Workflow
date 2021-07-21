def get_species_fasta(wildcards):
    return species_table.loc[wildcards.species]["pep_fasta"]


def get_longest_isoforms(wildcards):
    return os.path.join("FS/longest_isoforms/", os.path.split(get_species_fasta(wildcards))[1])


if config["auto_isoform_filtering"] == "YES":
    rule filter_isoforms:
        output:
            "FS/longest_isoforms/{species}.fa"
        params:
            fa = get_species_fasta,
            iso = get_longest_isoforms,
        shell:
            "python workflow/scripts/longest_isoforms.py {params.fa} && "
            "mv {params.iso} {output}"


    rule Orthofinder_link_all:
        input:
            "FS/longest_isoforms/{species}.fa",
        output:
            ORTHOFINDER + "{species}.fa"
        params: get_longest_isoforms,
        shell:
            "ln --symbolic $(readlink --canonicalize {input}) {output}"


else:
    rule Orthofinder_link_all_and_isoform_loc_copy:
        output:
            link = ORTHOFINDER + "{species}.fa",
            iso_link = "FS/longest_isoforms/{species}.fa"
        params:
            fa = get_species_fasta,
        shell:
            "ln --symbolic $(readlink --canonicalize {params.fa}) {output.link} && "
            "ln --symbolic $(readlink --canonicalize {params.fa}) {output.iso_link}"


rule Orthofinder_prepare:
    "Split fasta files, rename species and sequences and prepare blast databases"
    input:
        fastas = expand(ORTHOFINDER + "{species}.fa", species=SPECIES)
    output:
        check = touch(ORTHOFINDER + "prepare.check"),
        fastas = expand(ORTHOFINDER + "Species{species_number}.fa", species_number=[x for x in range(0, N_SPECIES)]),
        db = expand(ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd", database_number=[x for x in range(0, N_SPECIES)])
    params:
        fasta_dir = ORTHOFINDER,
        temp_dir1 = ORTHOFINDER + "OrthoFinder/Results_*",
        temp_dir2 = ORTHOFINDER + "OrthoFinder/Results_*/WorkingDirectory/"
    log:
        "logs/orthofinder/prepare.log"
    benchmark:
        "benchmarks/orthofinder/prepare.json"
#    conda:
#        "orthofinder.yml"
    shell:
        "orthofinder -t 16 "
         "--fasta {params.fasta_dir} "
         "--search diamond "
         "--only-prepare "
         "2> {log} 1>&2 && "
         "mv {params.temp_dir2}* {params.fasta_dir} && "
         "rm --recursive --force {params.temp_dir1}"


rule index_pep_FASTAs:
    input:
        ORTHOFINDER + "Species{species_number}.fa",
    output:
        ORTHOFINDER + "Species{species_number}.fa.fai",
    shell:
        "samtools faidx {input}"


# added choice of chunk usage or not ;D
with open('config/config.yaml', 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

if config["chunk_usage"] == "ON":
    rule Orthofinder_split:
        "Split the headers of the input pep fastas into multiple files"
        input:
            fai = ancient(ORTHOFINDER + "Species{species_number}.fa.fai")
        output:
            expand(ORTHOFINDER + "{{species_number}}/chunks/ids_{chunk_id}.tsv", chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)])
        params:
            folder = ORTHOFINDER,
            number_of_chunks = CHUNKS_ORTHO,
            species = "{species_number}"
        log:
        #send the log here as an exception; also controls the actual output 
            "orthofinder/{species_number}/split.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/split.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{params.folder}/{params.species}/chunks/ids_ "
            "2> {log}"


    rule Orthofinder_BLASTP:
        "Run blastp of each chunk"
        input:
            fasta = ORTHOFINDER + "Species{species_number}.fa",
            fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
            chunk = ORTHOFINDER + "{species_number}/chunks/ids_{chunk_id}.tsv",
            db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
        output:
            tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_{chunk_id}.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_{chunk_id}.json"
        threads: 6
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cut -f 1 {input.chunk} "
            "| xargs samtools faidx {input.fasta}"
            "| diamond blastp "
            "--ultra-sensitive "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 0.001 "
            "--out {output.tsv} "
            "--threads {threads} "
            "2> {log} 1>&2"



    rule Orthofinder_BLASTP_merge:
        "Merge results from the different blastps"
        input:
            expand(ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv", chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)])
        output:
            tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_merge.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_merge.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cat {input} > {output} 2> {log}"


## if chunks are not desired or feasible (anything other than "ON" under chunk_usage in the config.yaml)
else:
    rule Orthofinder_BLASTP:
        "Run blastp of each chunk"
        input:
            fasta = ORTHOFINDER + "Species{species_number}.fa",
            fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
            db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
        output:
            tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp.tsv"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp.json"
        threads: 6
    #    conda:
    #        "orthofinder.yml"
        shell:
            "diamond blastp --query {input.fasta} "
            "--ultra-sensitive "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 0.001 "
            "--out {output.tsv} "
            "--threads {threads} "
            "2> {log} 1>&2"


    rule Orthofinder_BLASTP_merge:
        "Merge results from the different blastps"
        input:
            expand(ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp.tsv")
        output:
            tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_merge.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_merge.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cat {input} > {output} 2> {log}"


rule Orthofinder_orthogroups:
    "Join blastp results, normalize bitscores and run mcl"
    input:
        tsv = expand(ORTHOFINDER + "Blast{database_number}_{species_number}.txt",
                  species_number = [x for x in range(0,N_SPECIES)],
                      database_number = [x for x in range(0,N_SPECIES)])
    output:
        touch(ORTHOFINDER + "groups.check"),
    params:
        fasta_dir = ORTHOFINDER,
        temp_dir = ORTHOFINDER + "OrthoFinder/Results_*",
#        inflation = params["orthofinder"]["mcl_inflation"],
    threads: 24 # There is no reason to go beyond this value - really?
    log:
        "logs/orthofinder/groups.log"
    benchmark:
        "benchmarks/orthofinder/groups.json"
#    conda:
#        "orthofinder.yml"
    shell:
        "orthofinder "
        "--algthreads {threads} "
        "--inflation 1.5 "
        "--blast {params.fasta_dir} "
        "--only-groups "
        "2> {log} 1>&2 && "
        "mv {params.temp_dir}/ {params.fasta_dir} "



rule Orthofinder_Trees:
    input:
        rules.Orthofinder_orthogroups.output
    output:
        touch(ORTHOFINDER + "trees.check")
    params:
        orthofinder_dir = ORTHOFINDER + "Results_*",
        tree_program = "fasttree",
        msa_program = "muscle",
#        mafft as alternative
#    conda:
#        "orthofinder.yml"
    log:
        "logs/orthofinder/trees.log"
    benchmark:
        "benchmarks/orthofinder/trees.bmk"
    threads: 24
    shell:
        "orthofinder "
        "--from-groups {params.orthofinder_dir} "
        "--only-trees "
        "--method msa "
        "--msa_program {params.msa_program} "
        "--tree_program fasttree "
        "--algthreads {threads} "
        "--threads {threads} "
        "2> {log} 1>&2"



rule Orthofinder_orthologues:
    input:
        ORTHOFINDER + "trees.check"
    output:
        touch(ORTHOFINDER + "orthologues.check")
#    conda: "orthofinder.yml"
    params:
        orthofinder_dir = ORTHOFINDER + "Results_*_1",
    log:
        "logs/orthofinder/orthologues.log"
    threads: 24
    shell:
        "orthofinder "
        "--from-trees {params.orthofinder_dir} "
        "-a {threads} "
        "--threads {threads} "
        "2> {log} 1>&2"



rule Orthofinder_cleanup:
    input:
        ORTHOFINDER + "orthologues.check"
    output:
        touch(ORTHOFINDER + "cleanup.check")
    params:
        orthofinder_dir = ORTHOFINDER,
        n_species = N_SPECIES - 1,
        results_dir_1 = "Results_*_1",
        results_dir_2 = "Results_*_2",
    log:
        "logs/orthofinder/clean.log"
    shell:
        "pushd {params.orthofinder_dir} 2> {log} 1>&2 && "
        "mkdir search/ && "
        "for i in {{0..{params.n_species}}} ; do "
        "mv ${{i}} search/ ; "
        "done && "
        "mv Blast* search/ &&"
        "mv diamond* search/ && "
        ""
        "mv {params.results_dir_1}/Gene_Trees/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/MultipleSequenceAlignments/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/Orthogroup_Sequences/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/WorkingDirectory/SpeciesTree_unrooted* {params.results_dir_2}/WorkingDirectory/ && "
        "mv {params.results_dir_1}/WorkingDirectory/*_ids/ {params.results_dir_2}/WorkingDirectory/ && "
        "mv {params.results_dir_2} final-results/ && "
        "rm -rf Results_*_1 && "
        "mv Results*/Comparative_Genomics_Statistics/* final-results/Comparative_Genomics_Statistics/ && "
        "mv Results*/Orthogroups/ final-results/ && "
        "rm -rf Results_* "


rule hypothesis_species_trees:
    input:
        ORTHOFINDER + "cleanup.check"
    output:
        "orthofinder/final-results/Species_Tree/hypothesis_specific/{hypothesis}/SpeciesTree_rooted_node_labels.txt",
    script:
        "../scripts/create_hypothesis_species_trees.R"


rule Orthofinder_complete:
    input:
        expand("orthofinder/final-results/Species_Tree/hypothesis_specific/{hypothesis}/SpeciesTree_rooted_node_labels.txt", hypothesis=HYPOTHESES)
#        ORTHOFINDER + "cleanup.check"
    output:
        touch(ORTHOFINDER + "complete.check")

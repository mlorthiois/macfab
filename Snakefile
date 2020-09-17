import glob

SOFTWARE=["stringtie", "flair"]
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all # never launch the all and the compact rules on the cluster

 # a simple rule to launch the full pipeline without specifiying the final rule (graph)
rule all:
    input:
        results="results/gffcompare/Graph.recap.pdf"
    threads:1
    resources: # is used by snakemake as input for the sbatch command
        ram="6G"


###############################################################################    
rule sam2bam:
    input:
        "results/minimap2/minimap.{annot}.sam"
    output:
        "results/minimap2/minimap.{annot}.sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads:10
    resources:
        ram="6G"
    shell:
        "samtools view -b -@ {threads} {input} | samtools sort -o {output}"


rule bam2bed12:
    input:
        "results/minimap2/minimap.{annot}.sorted.bam"
    output:
        "results/minimap2/minimap.{annot}.sorted.converted.bed12"
    conda:
        "envs/flair.yaml"
    shell:
        "bamToBed -bed12 -i {input} > {output}"


rule ungzip_genome_ref:
    input:
        config["reference_path"]
    output:
        "results/utilities/reference.dna.uncompressed.fa"
    run:
        print(input[0][-2:])
        if input[0][-2:]=="gz":
            shell("gunzip -c {input} > {output}")


rule gtfToBed12:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "results/gtfToBed12/{annot}.converted.bed12"
    conda:
        "envs/minimap.yaml" # defines a conda env to launch the command
    threads:1
    resources:
        ram="6G"
    shell:
        "paftools.js gff2bed {input} > {output}"


rule mapping:
    input:
        fastq=config["fastq_path"],
        bed="results/gtfToBed12/{annot}.converted.bed12",
        fa=config["reference_path"]
    output:
        "results/minimap2/minimap.{annot}.sam"
    conda:
        "envs/minimap.yaml"
    threads:10
    resources:
        ram="60G"
    shell:
        "minimap2 -t {threads} -ax splice --MD --junc-bed {input.bed} {input.fa} {input.fastq} > {output}"


###############################################################################
# a rule who installs a conda env with R, devtools et biocmanager
# the script installs inside the R env bambu 0.2 (release URL) and BSgenome (needed by bambu)
# ensure that the env is set before bambu, avoid conflicts when two bambu instances are launched at the same time
# rule configR:
#     output: touch("results/utilities/.R_config")
#     conda:
#         "envs/r.yaml"
#     script:
#         "scripts/install.R"


# analyze expression with bambu
# asks for .R_config created by the R_config rule
# a custom function is used to get the wildcards values ; see https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions
# rule bambu:
#     input:
#         gtf=lambda wildcards: config["annotation"][wildcards.annot],
#         bam="results/minimap2/minimap.{annot}.sorted.bam",
#         fa=config["reference_path"],
#         isConfig="results/utilities/.R_config"
#     output:
#         o_name="results/bambu/bambu.{annot}.gtf"
#     shadow: # snakemake symlinks the data, runs the soft inside the shadow dir then move the output out of it and delete it
#         "shallow"
#     conda:
#         "envs/r.yaml"
#     threads:1
#     resources:
#         ram="30G"
#     script:
#         "scripts/bambu.R"

###############################################################################
rule stringtie:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="results/minimap2/minimap.{annot}.sorted.bam"
    output:
        uncompressed_file="results/utilities/uncompress.{annot}.gtf",
        final_file="results/stringtie/stringtie.{annot}.gtf"
    conda:
        "envs/stringtie.yaml"
    threads:6
    resources:
        ram="10G"
    shell:
        """
        if [[ {input.gtf} =~ \.gz$ ]]; then 
            gunzip -c {input.gtf} > {output.uncompressed_file}
        fi  
        stringtie -L -G {output.uncompressed_file} -o {output.final_file} -p {threads} {input.bam}
         """



###############################################################################
rule flair_correct:
    input:
        query="results/minimap2/minimap.{annot}.sorted.converted.bed12",
        fa="results/utilities/reference.dna.uncompressed.fa",
        gtf="results/utilities/uncompress.{annot}.gtf"
    output:
        "results/flair/flair.{annot}_all_corrected.bed"
    params:
        "results/flair/flair.{annot}"
    conda:
        "envs/flair.yaml"
    shell:
        "flair.py correct -q {input.query} -g {input.fa} -f {input.gtf} -o {params}"


rule flair_collapse:
    input:
        fa="results/utilities/reference.dna.uncompressed.fa",
        reads=config["fastq_path"],
        query="results/flair/flair.{annot}_all_corrected.bed"
    output:
        "results/flair/flair.collapse.{annot}.isoforms.bed"
    params:
        "results/flair/flair.collapse.{annot}"
    conda:
        "envs/flair.yaml"
    shell:
        "flair.py collapse -g {input.fa} -r {input.reads} -q {input.query} -o {params}"


rule flairbedToGenePred:
    input:
        "results/flair/flair.collapse.{annot}.isoforms.bed"
    output:
        "results/flair/flair.{annot}.gtf"
    conda:
        "envs/bedToGenePred.yaml"
    params:
        "results/flair/flair.{annot}.gpf"
    shell:
        """
        bedToGenePred {input} {params} \
            | genePredToGtf file {params} {output}
        """




###############################################################################
# analyze expression using talon
# creates a dedicated config file on the fly using information provided in the config file   
# renames output to match the expected result  
# rule talon:
#     input:
#         fa=config["reference_path"],
#         sam="results/minimap2/minimap.{annot}.sam",
#         gtf=lambda wildcards: config["annotation"][wildcards.annot]
#     output:
#         gtf="results/talon/talon.{annot}.gtf",
#         db="results/talon/talon.{annot}.db"
#     conda:
#         "envs/talon.yaml"
#     threads:20
#     priority: 50 # ensure that the rule is created top priority to avoid disk space issues
#     shadow: "shallow"
#     resources:
#         ram="50G"
#     params:
#         cell_line=config["cell_line"],
#         prefix="talon.{annot}",
#         used_annot="{annot}",
#         reference_build=config["reference_build"],
#         data_type=config["data_type"]
#     shell:
#         """
#         ~/.local/bin/talon_label_reads --f {input.sam} --g {input.fa} --o {params.prefix} --t={threads}
#         ~/.local/bin/talon_initialize_database --f {input.gtf} --g {params.reference_build} --a {params.used_annot} --idprefix {params.prefix} --o {params.prefix}
#         echo {params.cell_line},{params.data_type},nanopore,{params.prefix}_labeled.sam > talon.config
#         ~/.local/bin/talon --f talon.config --db {output.db} --build {params.reference_build} -t {threads} --o {threads}
#         ~/.local/bin/talon_create_GTF --db {output.db} -b {params.reference_build} -a {params.used_annot} --o {output.gtf}
#         mv {output.gtf}_talon.gtf {output.gtf}
#         """


###############################################################################
# intersect gtf with bam to keep only expressed genes and not those from the annotation        
rule only_seen_exons:
    input:
        gtf="results/{software}/{software}.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam"
    output:
        final="results/gffcompare/{software}.{annot}.filtered.gtf",
        exon=temp("{software}.{annot}.EO.gtf")
    threads:1
    conda:
        "envs/bedtools.yaml"
    resources:
        ram="50G"
    shell:
        """
        grep $'\t'exon$'\t' {input.gtf} > {output.exon}
        bedtools intersect -u -s -split -a {output.exon} -b {input.bam} > {output.final}
        """

# compare annotation to ref      
rule gffcompare:
    input:
        test="results/gffcompare/{software}.{annot}.filtered.gtf",
        ref="results/utilities/uncompress.{annot}.gtf"
    output:
        result="results/gffcompare/{software}.{annot}.stats"
    threads:1
    conda:
        "envs/gffcompare.yaml"
    shadow: "shallow"
    params:
        prefix="{software}_{annot}",
    resources:
        ram="6G"
    shell:
         """
         gffcompare {input.test} -r {input.ref} -o {params.prefix}
         mv {params.prefix}.stats {output.result}
         """


# format values to use it as R input      
rule parse_gffcompare:
    input:
        expand("results/gffcompare/{software}.{annot}.stats", software=SOFTWARE, annot=config["annotation"])
    output:
        Sensitivity="results/gffcompare/Sensitivity.gffparse.tsv",
        Values="results/gffcompare/Values.gffparse.tsv"
    threads:1
    resources:
        ram="6G"
    script:
        "scripts/gffparse.py"


# plot the expected results with ggplot2
rule graph:
    input:
        Sensitivity="results/gffcompare/Sensitivity.gffparse.tsv",
        Values="results/gffcompare/Values.gffparse.tsv"
    output:
        "results/gffcompare/Graph.recap.pdf"
    threads:1
    conda:
        "envs/r.yaml"
    resources:
        ram="6G"
    script:
        "scripts/graphs.R"

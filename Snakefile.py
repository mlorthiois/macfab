import glob

# rule mapping
SOFTWARE=["stringtie", "talon", "flair", "bambu"] # NEVER add annotating rules without adding software here
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all, compact, configR # never launch the all and the compact rules on the cluster

 # a simple rule to launch the full pipeline without specifiying the final rule (graph)
rule all:
    input:
        results="Graph.recap.pdf"
    threads:1
    resources: # is used by snakemake as input for the sbatch command
        ram="6G"

# a rule who installs a conda env with R, devtools et biocmanager
# the script installs inside the R env bambu 0.2 (release URL) and BSgenome (needed by bambu)
# ensure that the env is set before bambu, avoid conflicts when two bambu instances are launched at the same time
rule configR:
    output: touch(".R_config")
    conda:
        "envs/r.yaml"
    script:
        "scripts/install.R"

# a rule who creates one fastq file from the dir
# bambu doesn't know how to use a dir
rule compact:
    output:
       temp("compacted.fastq")
    params:
        config["fastq_path"]
    shell:
        "cat {params} > {output}"

# convert annotations for minimap
# custom path because at the time of the writing the conda paftools.js doesn't work
rule gtfToBed12:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        temp("{annot}.converted.bed12")
    conda:
        "envs/minimap.yaml" # defines a conda env to launch the command
    threads:1
    resources:
        ram="6G"
    shell:
        config["paftools.js"] + " gff2bed {input} > {output}"

# map fastq on the ref genome with annotations provided
rule mapping:
    input:
        fastq=rules.compact.output,
        bed="{annot}.converted.bed12",
        fa=config["reference_path"]
    output:
        temp("minimap.{annot}.sam")
    conda:
        "envs/minimap.yaml"
    threads:10
    resources:
        ram="60G"
    shell:
        "minimap2 -t {threads} -ax splice --MD --junc-bed {input.bed} {input.fa} {input.fastq} > {output}"

# convert sam to bam and sort it      
rule sam2bam:
    input:
        "minimap.{annot}.sam"
    output:
        "minimap.{annot}.sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads:10
    resources:
        ram="6G"
    shell:
        "samtools view -b -@ {threads} {input} | samtools sort -o {output}"

# analyze expression with bambu
# asks for .R_config created by the R_config rule
# a custom function is used to get the wildcards values ; see https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions
rule bambu:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap.{annot}.sorted.bam",
        fa=config["reference_path"],
        isConfig=".R_config"
    output:
        o_name="bambu.{annot}.gtf"
    shadow: # snakemake symlinks the data, runs the soft inside the shadow dir then move the output out of it and delete it
        "shallow"
    conda:
        "envs/r.yaml"
    threads:1
    resources:
        ram="30G"
    script:
        "scripts/bambu.R"

# analyze expression using stringtie
# stringtie needs to be installed and the path specified using the config file
# no shadow dir here because stringtie only creates its output
rule stringtie:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap.{annot}.sorted.bam"
    output:
        "stringtie.{annot}.gtf"
    threads:6
    resources:
        ram="10G"
    shell:
        config["stringtie"] + " -L -G {input.gtf} -o {output} -p {threads} {input.bam}"

# annotates using flair
# flair needs bed12 mapping
# flair outputs bed who needs to be converted to gtf
rule flair:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap.{annot}.sorted.bam",
        fa=config["reference_path"],
        fastq="compacted.fastq"
    output:
        "flair.{annot}.gtf"
    conda:
        "envs/flair.yaml"
    shadow: "shallow"
    threads:40
    resources:
        ram="30G"
    params:
        bed12=config["bed12ToGtf"], # needed to use the value in multi-lines shell command
        prefix="flair.{annot}" # software use prefix but prefix can't be used as output (because no file matching exactly this name will be created)
    shell:
        """
        python3 -m pip install kerneltree
        python3 -m pip install Cython
        bamToBed -bed12 -i {input.bam} > converted.bed12
        flair.py correct -q converted.bed12 -g {input.fa} -f {input.gtf} -o {params.prefix}
        flair.py collapse -g {input.fa} -r {input.fastq} -q {params.prefix}_all_corrected.bed -o {params.prefix}
        {params.bed12} {params.prefix}.isoforms.bed > {output}
        """
# analyze expression using talon
# creates a dedicated config file on the fly using information provided in the config file   
# renames output to match the expected result  
rule talon:
    input:
        fa=config["reference_path"],
        sam="minimap.{annot}.sam",
        gtf=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        gtf="talon.{annot}.gtf",
        db="talon.{annot}.db"
    conda:
        "envs/talon.yaml"
    threads:20
    priority: 50 # ensure that the rule is created top priority to avoid disk space issues
    shadow: "shallow"
    resources:
        ram="50G"
    params:
        cell_line=config["cell_line"],
        prefix="talon.{annot}",
        used_annot="{annot}",
        reference_build=config["reference_build"],
        data_type=config["data_type"]
    shell:
        """
        ~/.local/bin/talon_label_reads --f {input.sam} --g {input.fa} --o {params.prefix} --t={threads}
        ~/.local/bin/talon_initialize_database --f {input.gtf} --g {params.reference_build} --a {params.used_annot} --idprefix {params.prefix} --o {params.prefix}
        echo {params.cell_line},{params.data_type},nanopore,{params.prefix}_labeled.sam > talon.config
        ~/.local/bin/talon --f talon.config --db {output.db} --build {params.reference_build} -t {threads} --o {threads}
        ~/.local/bin/talon_create_GTF --db {output.db} -b {params.reference_build} -a {params.used_annot} --o {output.gtf}
        mv {output.gtf}_talon.gtf {output.gtf}
        """

# intersect gtf with bam to keep only expressed genes and not those from the annotation        
rule only_seen_exons:
    input:
        gtf="{software}.{annot}.gtf",
        bam="minimap.{annot}.sorted.bam"
    output:
        final="{software}.{annot}.filtered.gtf",
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
        test="{software}.{annot}.filtered.gtf",
        ref=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        result="{software}.{annot}.stats"
    threads:1
    shadow: "shallow"
    params:
        prefix="{software}_{annot}",
        gffcompare=config["gffcompare"]
    resources:
        ram="6G"
    shell:
         """
         {params.gffcompare} {input.test} -r {input.ref} -o {params.prefix}
         mv {params.prefix}.stats {output.result}
         """
# format values to use it as R input      
rule parse_gffcompare:
    input:
        expand("{software}.{annot}.stats", software=SOFTWARE, annot=config["annotation"].keys())
    output:
        Sensitivity="Sensitivity.gffparse.tsv",
        Values="Values.gffparse.tsv"
    threads:1
    resources:
        ram="6G"
    script:
        "scripts/gffparse.py"

# plot the expected results with ggplot2
rule graph:
    input:
        Sensitivity="Sensitivity.gffparse.tsv",
        Values="Values.gffparse.tsv"
    output:
        "Graph.recap.pdf"
    threads:1
    conda:
        "envs/r.yaml"
    resources:
        ram="6G"
    script:
        "scripts/graphs.R"

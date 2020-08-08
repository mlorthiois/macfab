import glob

# rule mapping
SOFTWARE=["stringtie", "talon", "flair", "bambu"] # NEVER add annotating rules without adding software here
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all, compact, configR # never launch the all and the compact rules on the cluster

rule all: # a simple rule to launch the full pipeline without specifiying the final rule (graph)
    input:
        results="Graph.recap.pdf"
    threads:1
    resources: # is used by snakemake as input for the sbatch command
        ram="6G"

rule configR:
    output: touch(".R_config")
    conda:
        "envs/r.yaml"
    script:
        "scripts/install.R"

rule compact:
    input:
        lambda files: glob.glob(config["fastq_path"])
    output:
       "compacted.fastq"
    threads:1
    resources:
        ram="6G"
    shell:
        "cat {input} > {output}"

rule gtfToBed12:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "{annot}.converted.bed12"
    conda:
        "envs/minimap.yaml" # defines a conda env to launch the command
    threads:1
    resources:
        ram="6G"
    shell:
        config["paftools.js"] + " gff2bed {input} > {output}"

rule mapping:
    input:
        fastq=rules.compact.output,
        bed="{annot}.converted.bed12",
        fa=config["reference_path"]
    output:
        "minimap.{annot}.sam"
    conda:
        "envs/minimap.yaml"
    threads:10
    resources:
        ram="20G"
    shell:
        "minimap2 -t {threads} -ax splice --MD --junc-bed {input.bed} {input.fa} {input.fastq} > {output}"
        
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

rule bambu:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap.{annot}.sorted.bam",
        fa=config["reference_path"],
        isConfig=".R_config"
    output:
        o_name="bambu.{annot}.gtf"
    shadow:
        "shallow"
    conda:
        "envs/r.yaml"
    threads:1
    resources:
        ram="10G"
    script:
        "scripts/bambu.R"
        
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
    threads:20
    resources:
        ram="10G"
    params:
        bed12=config["bed12ToGtf"], # needed to use the value in multi-lines shell command
        prefix="flair.{annot}" # software use prefix but prefix can't be used as output (because no file matching exactly this name will be created)
    shell:
        """
        bamToBed -bed12 -i {input.bam} > converted.bed12
        flair.py correct -q converted.bed12 -g {input.fa} -f {input.gtf} -o {params.prefix}
        flair.py collapse -g {input.fa} -r {input.fastq} -q {params.prefix}_all_corrected.bed -o {params.prefix}
        {params.bed12} {params.prefix}.isoforms.bed > {output}
        """
        
rule talon:
    input:
        fa=config["reference_path"],
        sam="minimap.{annot}.sam",
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        isConfig=".shell_config"
    output:
        gtf="talon.{annot}.gtf",
        db="talon.{annot}.db"
    conda:
        "envs/talon.yaml"
    threads:20
    shadow: "shallow"
    resources:
        ram="50G"
    params:
        cell_line=config["cell_line"],
        prefix="talon.{annot}",
        used_annot="{annot}"
    shell:
        """
        ~/.local/bin/talon_label_reads --f {input.sam} --g {input.fa} --o {params.prefix} --t={threads}
        ~/.local/bin/talon_initialize_database --f {input.gtf} --g CanFam3 --a {params.used_annot} --idprefix {params.prefix} --o {params.prefix}
        echo {params.cell_line},Dog_transcript,nanopore,{params.prefix}_labeled.sam > talon.config
        ~/.local/bin/talon --f talon.config --db {output.db} --build CanFam3 -t {threads} --o {threads}
        ~/.local/bin/talon_create_GTF --db {output.db} -b CanFam3 -a {params.used_annot} --o {output.gtf}
        mv {output.gtf}_talon.gtf {output.gtf}
        """
        
rule only_seen_exons:
    input:
        gtf="{software}.{annot}.gtf",
        fastq="compacted.fastq",
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
        bedtools intersect -wa -s -split -a {output.exon} -b {input.bam} > {output.final}
        """
        
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
 
rule graph:
    input:
        Sensitivity="Sensitivity.gffparse.tsv",
        Values="Values.gffparse.tsv"
    output:
        "Graph.recap.pdf"
    threads:1
    resources:
        ram="6G"
    script:
        "scripts/graphs.R"

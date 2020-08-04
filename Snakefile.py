import glob

# rule mapping
SOFTWARE=["stringtie", "talon", "flair", "bambu"] # NEVER add annotating rules without adding software here
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all, compact # never launch the all and the compact rules on the cluster

rule all: # a simple rule to launch the full pipeline without specifiying the final rule (graph)
    input:
        results="Graph.recap.pdf",
        config_done=".config_done"
    threads:1
    resources: # is used by snakemake as input for the sbatch command
        ram="6G"

rule config:
    output: touch(".config_done")
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
    log: "gtf_{annot}.log"
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
    log: "mapping_{annot}.log"
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
    log: "sam2bam_{annot}.log"
    resources:
        ram="6G"
    shell:
        "samtools view -b -@ {threads} {input} | samtools sort -o {output}"

rule bambu:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap.{annot}.sorted.bam",
        fa=config["reference_path"]
    output:
        o_dir=directory("bambu.{annot}"),
        o_name="bambu.{annot}.gtf"
    conda:
        "envs/r.yaml"
    threads:1
    log: "bambu_{annot}.log"
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
    log: "stringtie_{annot}.log"
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
    threads:1
    resources:
        ram="10G"
    log: "flair_{annot}.log"
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
        gtf=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "talon.{annot}.gtf"
    conda:
        "envs/talon.yaml"
    threads:10
    resources:
        ram="20G"
    log: "talon_{annot}.log"
    params:
        cell_line=config["cell_line"],
        prefix="talon.{annot}",
        used_annot="{annot}"
    shell:
        #FIXME ensure that the database isn't created twice
        """
        talon_label_reads --deleteTmp --f {input.sam} --g {input.fa} --o {params.prefix} --t={threads}
        if [ -f "{params.prefix}.db" ]; then
        rm "{params.prefix}.db"
        fi
        talon_initialize_database --f {input.gtf} --g CanFam3 --a {params.used_annot} --idprefix {params.prefix} --o {params.prefix}
        echo {params.cell_line},Dog_transcript,nanopore,{params.prefix}_labeled.sam > talon.config
        talon --f talon.config --db {params.prefix}.db --build CanFam3 -t {threads} --o {threads}
        talon_create_GTF --db {params.prefix}.db -b CanFam3 -a {params.used_annot} --o {output}
        """
        
rule only_seen_exons:
    input:
        gtf="{software}.{annot}.gtf",
        fastq="compacted.fastq",
        bam="minimap.{annot}.sorted.bam"
    output:
        final="{software}.{annot}.filtered.gtf",
        exon=temp("{software}.{annot}.EO.gtf"),
        sorted_gtf=temp("{software}.{annot}.EO.sorted.gtf")
    threads:1
    log: "only_seen_exons_{annot}_{software}.log"
    conda:
        "envs/bedtools.yaml"
    resources:
        ram="20G"
    shell:
        """
        grep $'\t'exon$'\t' {input} > {output.exon}
        sort -k1,1 -k4,4 {output.exon} > {output.sorted_gtf}
        bedtools intersect -sorted -wa -s -split -a {output.sorted_gtf} -b {input.bam} > {output.final}
        """
        
rule gffcompare:
    input:
        test="{software}.{annot}.filtered.gtf",
        ref=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        folder=directory("{software}.{annot}.folder"),
        result="{software}.{annot}.stats"
    threads:1
    log: "gffcompare_{annot}_{software}.log"
    resources:
        ram="6G"
    shell:
        config["gffcompare"] + " {input.test} -r {input.ref} -o {output.folder}"
        
rule parse_gffcompare:
    input:
        expand("{software}.{annot}.stats", software=SOFTWARE, annot=config["annotation"].keys())
    output:
        Sensitivity="Sensitivity.parsed.tsv",
        Values="Values.parsed.tsv"
    threads:1
    log: "parse_gffcompare.log"
    resources:
        ram="6G"
    script:
        "scripts/gfffparse.py"
 
rule graph:
    input:
        Sensitivity="Sensitivity.parsed.tsv",
        Values="Values.parsed.tsv"
    output:
        "Graph.recap.pdf"
    log: "graph.log"
    threads:1
    resources:
        ram="6G"
    script:
        "scripts/graph.R"

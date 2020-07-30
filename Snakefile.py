# rule mapping
SOFTWARE=["stringtie", "talon", "flair", "bambu"]

configfile: "config.yaml"
localrules: all, compact

rule all:
    input:
        "Sensitivity_parsed.tsv"
    threads:1
    params:
        ram="6G"

rule compact:
    input:
        config["fastq_path"]
    output:
       "compacted.fastq"
    threads:1
    params:
        ram="6G"
    shell:
        "cat {input} > {output}"

rule gtfToBed12:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "{annot}_converted.bed12"
    conda:
        "envs/minimap.yaml"
    threads:1
    params:
        ram="6G"
    shell:
        "paftools.js gff2bed {input}"

rule mapping:
    input:
        fastq=rules.compact.output,
        bed="{annot}_converted.bed12",
        fa=config["reference_path"]
    output:
        "minimap_{annot}.sam"
    conda:
        "envs/minimap.yaml"
    threads:10
    params:
        ram="20G"
    shell:
        "minimap2 -t {threads} -ax splice --MD --junc-bed {input.bed} {intput.fa} {input.fastq} > {output}"
        
rule bam2sam:
    input:
        "minimap_{annot}.sam"
    output:
        "minimap_{annot}_sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads:10
    params:
        ram="6G"
    shell:
        "samtools view -b -@ {threads} {input} | samtools sort -o {output}"

rule bambu:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam=rules.bam2sam.output,
        fa=config["reference_path"]
    output:
        o_dir=directory("bambu_annot/"),
        o_name="bambu_{annot}.gtf"
    threads:1
    params:
        ram="10G"
    script:
        "scripts/bambu.R"
        
rule stringtie:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap_{annot}_sorted.bam"
    output:
        "stringtie_{annot}.gtf"
    threads:6
    params:
        ram="10G"
    shell:
        config["stringtie"] + " -L -G {input.gtf} -o {output} -p {threads} {input.bam}"

rule flair:
    input:
        gtf=lambda wildcards: config["annotation"][wildcards.annot],
        bam="minimap_{annot}_sorted.bam",
        fa=config["reference_path"],
        fastq="compacted.fastq"
    output:
        o_prefix="flair_{reference_path}",
        o_gtf="flair_{annot}.gtf"
    conda:
        "envs/flair.yaml"
    threads:1
    params:
        ram="10G"
    run:
        shell("bamToBed -bed12 -i {input.bam} > converted.bed12")
        shell("flair.py correct -q converted.bed12 -g {input.fa} -f {input.gtf} -o {output.o_prefix}")
        shell("flair.py collapse -g {input.fa} -r {input.fastq} -q {output.o_prefix}_all_corrected.bed -o {output.o_prefix}")
        shell(config["bed12togtf"] + " {output.o_prefix}.collapse.isoforms.bed > {output.o_gtf}")
        
rule talon:
    input:
        fa=config["reference_path"],
        sam="minimap_{annot}.sam",
        gtf=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        o_gtf="talon_{annot}.gtf",
        o_prefix="talon_{annot}"
    threads:10
    params:
        ram="20G"
    run:
        shell("talon_label_reads --f {input.sam} --g {input.fa} --o {output.o_prefix} --t={threads}")
        shell("talon_initialize_database --f {input.gtf} --g CanFam3 --a {annot} --idprefix {o.prefix} --o {output.o_prefix}")
        shell("cat " + config["cell_line"] + ",Dog_transcript,nanopore,{output.o_prefix}_labelled.sam > talon.config")
        shell("talon --f talon.config --db {output.o_prefix}.db --build CanFam3 -t {threads} --o {threads}")
        shell("talon_create_GTF --db {output.o_prefix}.db -b CanFam3 -a {annot} --o {output.o_prefix}")

#FIXME add gtf intersect/filtration time
        
rule gffcompare:
    input:
        test="{software}_{annot}.gtf",
        ref=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "{software}_{annot}"
    threads:1
    params:
        ram="6G"
    shell:
        config["gffcompare"] + " {input.test} -r {input.ref} -o {output}"
        
rule parse_gffcompare:
    input:
        expand("{software}_{annot}.stats", software=SOFTWARE, annot=config["annotation"].keys())
    output:
        Sensitivity="Sensitivity_parsed.tsv",
        Values="Values_parsed.tsv"
    threads:1
    params:
        ram="6G"
    script:
        "scripts/gfffparse.py"
 
rule graph:
    input:
        Sensitivity="Sensitivity_parsed.tsv",
        Values="Values_parsed.tsv"
    output:
        "Graph_recap.pdf"
    threads:1
    params:
        ram="6G"
    script:
        "scripts/graph.R"

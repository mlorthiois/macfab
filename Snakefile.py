import glob
# rule mapping
SOFTWARE=["stringtie", "talon", "flair", "bambu"]
wildcard_constraints:
    annot="[^./]+"

configfile: "config.yaml"
localrules: all, compact

rule all:
    input:
        "Graph.recap.pdf"
    threads:1
    resources:
        ram="6G"

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
        "envs/minimap.yaml"
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
        fa=config["reference_path"]
    output:
        o_dir=directory("bambu.{annot}"),
        o_name="bambu.{annot}.gtf"
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
        o_prefix="flair.{annot}",
        o_gtf="flair.{annot}.gtf"
    conda:
        "envs/flair.yaml"
    threads:1
    resources:
        ram="10G"
    shell:
        """
        bamToBed -bed12 -i {input.bam} > converted.bed12
        flair.py correct -q converted.bed12 -g {input.fa} -f {input.gtf} -o {output.o_prefix}
        flair.py collapse -g {input.fa} -r {input.fastq} -q {output.o_prefix}_all_corrected.bed -o {output.o_prefix}
        ~tderrien/bin/convert/bed12Togtf.sh {output.o_prefix}.collapse.isoforms.bed > {output.o_gtf}
        """
        
rule talon:
    input:
        fa=config["reference_path"],
        sam="minimap.{annot}.sam",
        gtf=lambda wildcards: config["annotation"][wildcards.annot]
    output:
        o_gtf="talon.{annot}.gtf",
        o_prefix="talon.{annot}"
    threads:10
    resources:
        ram="20G"
    run:
        annot=wildcards.annot
        shell("talon_label_reads --f {input.sam} --g {input.fa} --o {output.o_prefix} --t={threads}")
        shell("talon_initialize_database --f {input.gtf} --g CanFam3 --a " + annot + " --idprefix {o.prefix} --o {output.o_prefix}")
        shell("cat " + config["cell_line"] + ",Dog_transcript,nanopore,{output.o_prefix}_labelled.sam > talon.config")
        shell("talon --f talon.config --db {output.o_prefix}.db --build CanFam3 -t {threads} --o {threads}")
        shell("talon_create_GTF --db {output.o_prefix}.db -b CanFam3 -a " + annot + " --o {output.o_prefix}")

rule only_seen_exons:
    input:
        gtf="{software}.{annot}.gtf",
        fastq="compacted.fastq"
        
    output:
        final="{software}.{annot}.filtered.gtf",
        exon=temp("{software}.{annot}.EO.gtf"),
        sorted_gtf=temp("{software}.{annot}.EO.sorted.gtf")
    threads:1
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
        folder="{software}.{annot}.folder",
        result="{software}.{annot}.stats"
    threads:1
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
    threads:1
    resources:
        ram="6G"
    script:
        "scripts/graph.R"

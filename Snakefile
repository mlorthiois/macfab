import glob

SOFTWARE=["stringtie", "talon", "flair", "bambu"]
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all # never launch the all and the compact rules on the cluster

 # a simple rule to launch the full pipeline without specifiying the final rule (graph)
rule all:
    input:
        gffcompare_recap="results/gffcompare/Graph.recap.pdf",
        sqanti_summary = "results/SQANTI3/SQANTI_report.pdf",
        rseqc = expand("results/RSeQC/{annot}.geneBodyCoverage.curves.pdf", annot=config["annotation"]),
        bai = expand("results/minimap2/minimap.{annot}.sorted.bam.bai", annot=config["annotation"])


###############################################################################
rule sam2bam:
    input:
        "results/minimap2/minimap.{annot}.sam"
    output:
        "results/minimap2/minimap.{annot}.sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads:30
    resources:
        ram="16G"
    log: "logs/{annot}_sam2bam.log"
    shell:
        """
        samtools view -b -@ {threads} {input} | samtools sort -o {output}
        """


rule bam2bai:
    input:
        "results/minimap2/minimap.{annot}.sorted.bam"
    output:
        "results/minimap2/minimap.{annot}.sorted.bam.bai"
    conda:
        "envs/samtools.yaml"
    threads:16
    resources:
        ram="16G"
    log: "logs/{annot}_sam2bam.log"
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule bam2bed12:
    input:
        "results/minimap2/minimap.{annot}.sorted.bam"
    output:
        temp("results/minimap2/minimap.{annot}.sorted.converted.bed12")
    conda:
        "envs/flair.yaml"
    log: "logs/{annot}_bam2bed12.log"
    threads:2
    resources:
        ram="16G"
    shell:
        "bamToBed -bed12 -i {input} > {output} 2> {log}"


rule ungzip_genome_ref:
    input:
        config["reference_path"]
    output:
        temp("results/utilities/reference.dna.uncompressed.fa")
    log: "logs/unzip_genome_ref.log"
    threads:1
    resources:
        ram="10G"
    run:
        if input[0][-2:]=="gz":
            shell("gunzip -c {input} > {output}")
        else:
            shell("cp {input} {output}")


rule ungzip_query_fastq:
    input:
        config["fastq_path"]
    output:
        temp("results/utilities/Query.uncompressed.fastq")
    log: "logs/unzip_query_fastq.log"
    threads:1
    resources:
        ram="16G"
    run:
        if input[0][-2:]=="gz":
            shell("gunzip -c {input} > {output}")
        else:
            shell("cp {input} {output}")


rule ungzip_gtf:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        temp("results/utilities/uncompress.{annot}.gtf")
    log: "logs/{annot}_unzip_gtf.log"
    threads:1
    resources:
        ram="10G"
    run:
        if input[0][-2:]=="gz":
            shell("gunzip -c {input} > {output}")
        else:
            shell("cp {input} {output}")


rule install_paftools:
    # Use custom paftools because last release doesn't work with GenBank GTF
    threads:1
    resources:
        ram="2G"
    output:
        "results/utilities/paftools.js"
    log: "logs/install_paftools.log"
    shell:
        """
        curl -L https://raw.githubusercontent.com/lh3/minimap2/58c2251b18e70cdaa2e8e2088899001cfe7d69ae/misc/paftools.js -o {output} 2> {log}
        chmod +x {output} 2> {log}
        """


rule gtfToBed12:
    input:
        paftools = "results/utilities/paftools.js",
        gtf = lambda wildcards: config["annotation"][wildcards.annot]
    output:
        temp("results/utilities/{annot}.converted.bed12")
    conda:
        "envs/minimap.yaml"
    log: "logs/{annot}_gtfToBed12.log"
    threads:1
    resources:
        ram="8G"
    shell:
        """
        {input.paftools} gff2bed {input.gtf} > {output} 2> {log}
        """


###############################################################################
rule mapping:
    input:
        fastq="results/utilities/Query.uncompressed.fastq",
        bed="results/utilities/{annot}.converted.bed12",
        fa=config["reference_path"]
    output:
        temp("results/minimap2/minimap.{annot}.sam")
    conda:
        "envs/minimap.yaml"
    log: "logs/{annot}_mapping.log"
    threads: 30
    resources:
        ram="40G"
    shell:
        """
        minimap2 -t {threads} -ax splice --MD \
        --junc-bed {input.bed} {input.fa} {input.fastq} > {output} 2> {log}
        """


###############################################################################
rule bambu:
    input:
        gtf="results/utilities/uncompress.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam",
        fa="results/utilities/reference.dna.uncompressed.fa"
    output:
        bambu_gtf="results/bambu/bambu.{annot}.gtf"
    shadow: # snakemake symlinks the data, runs the soft inside the shadow dir then move the output out of it and delete it
        "shallow"
    log: "logs/{annot}_bambu.log"
    conda:
        "envs/r.yaml"
    threads:2
    resources:
        ram="24G"
    script:
        "scripts/bambu.R"


###############################################################################
rule stringtie:
    input:
        gtf="results/utilities/uncompress.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam"
    output:
        final_file="results/stringtie/stringtie.{annot}.gtf"
    conda:
        "envs/stringtie.yaml"
    threads: 16
    resources:
        ram="20G"
    log: "logs/{annot}_stringtie.log"
    shadow:
        "shallow"
    shell:
        """
        stringtie -L -G {input.gtf} -o {output.final_file} -p {threads} {input.bam} 2> {log}
        """


###############################################################################
rule flair_correct:
    input:
        query="results/minimap2/minimap.{annot}.sorted.converted.bed12",
        fa="results/utilities/reference.dna.uncompressed.fa",
        gtf="results/utilities/uncompress.{annot}.gtf"
    output:
        temp("results/flair/flair.{annot}_all_corrected.bed")
    params:
        "results/flair/flair.{annot}"
    log: "logs/{annot}_flair_correct.log"
    threads:2
    resources:
        ram="16G"
    conda:
        "envs/flair.yaml"
    shell:
        "flair.py correct -q {input.query} -g {input.fa} -f {input.gtf} -o {params} 2> {log}"


rule flair_collapse:
    input:
        fa="results/utilities/reference.dna.uncompressed.fa",
        reads="results/utilities/Query.uncompressed.fastq",
        query="results/flair/flair.{annot}_all_corrected.bed",
        ref_gtf = "results/utilities/uncompress.{annot}.gtf"
    output:
        "results/flair/flair.{annot}.gtf"
    params:
        "results/flair/flair.collapse.{annot}"
    log:
        "logs/{annot}_flair_collapse.log"
    threads: 16
    resources:
        ram="20G"
    conda:
        "envs/flair.yaml"
    shell:
        """
        flair.py collapse -t {threads} -g {input.fa} -r {input.reads} -q {input.query} -o {params} -f {input.ref_gtf}
        mv {params}.isoforms.gtf {output} 2> {log}
        """


###############################################################################
rule talon_initialize_database:
    input:
        gtf="results/utilities/uncompress.{annot}.gtf"
    output:
        db="results/talon/talon.{annot}.db"
    shadow:
        "shallow"
    params:
        reference_genome_name = config["reference_genome_name"],
        prefix = "results/talon/talon.{annot}",
        annotation_name = "{annot}_annot"
    threads:2
    resources:
        ram="10G"
    conda:
        "envs/talon.yaml"
    log: "logs/{annot}_talon_initialize_database.log"
    shell:
        """
        talon_initialize_database --f {input.gtf} \
            --g {params.reference_genome_name} \
            --a {params.annotation_name} \
            --o {params.prefix} 2> {log}
        """


rule talon_label_reads:
    input:
        sam="results/minimap2/minimap.{annot}.sam",
        fa="results/utilities/reference.dna.uncompressed.fa"
    output:
        temp("results/talon/talon.{annot}_labeled.sam")
    params:
        prefix="results/talon/talon.{annot}"
    conda:
        "envs/talon.yaml"
    log: "logs/{annot}_talon_label_reads.log"
    shadow:
        "shallow"
    threads: 16
    resources:
        ram="30G"
    shell:
        """
        talon_label_reads --f {input.sam} \
            --g {input.fa} \
            --t {threads} \
            --o {params.prefix} 2> {log}
        """


rule create_talon_configfile:
    threads:1
    resources:
        ram="2G"
    input:
    output:
        "results/talon/talon.{annot}.config.csv"
    params:
        id=config["sample_id"],
        description=config["sample_description"],
        ngs="Nanoseq",
        sam="results/talon/talon.{annot}_labeled.sam"
    shadow:
        "shallow"
    log: "logs/{annot}_create_talon_configfile.log"
    shell:
       "echo {params.id},{params.description},{params.ngs},{params.sam} > {output} 2> {log}"


rule talon:
    input:
        config="results/talon/talon.{annot}.config.csv",
        db="results/talon/talon.{annot}.db",
        sam="results/talon/talon.{annot}_labeled.sam"
    output:
        "results/talon/talon.{annot}_QC.log"
    params:
        build=config["reference_genome_name"],
        prefix="results/talon/talon.{annot}"
    conda:
        "envs/talon.yaml"
    log: "logs/{annot}_talon.log"
    shadow:
        "shallow"
    threads:30
    resources:
        ram="50G"
    shell:
        """
        talon --f {input.config} \
            --db {input.db} \
            --build {params.build} \
            --threads {threads} \
            --o {params.prefix} 2> {log}
        """


rule filter_transcripts:
    input:
        db="results/talon/talon.{annot}.db",
        log="results/talon/talon.{annot}_QC.log",
        annot="results/utilities/uncompress.{annot}.gtf"
    output:
        "results/talon/talon.{annot}_filtered_transcripts.csv"
    log: "logs/talon_filter_transcripts_{annot}.log"
    params:
        annot="{annot}_annot",
        dataset=config["sample_id"]
    shadow:
        "shallow"
    conda:
        "envs/talon.yaml"
    threads:1
    resources:
        ram="20G"
    shell:
        """
        talon_filter_transcripts \
            --db {input.db} \
            -a {params.annot} \
            --o {output} 2> {log}
        """


rule talon_create_GTF:
    input:
        db="results/talon/talon.{annot}.db",
        csv="results/talon/talon.{annot}_filtered_transcripts.csv"
    output:
        "results/talon/talon.{annot}.gtf"
    params:
        annotation_name="{annot}_annot",
        ref_name=config["reference_genome_name"],
        prefix="results/talon/talon.{annot}"
    log: "logs/{annot}_talon_create_GTF.log"
    shadow:
        "shallow"
    conda:
        "envs/talon.yaml"
    threads:1
    resources:
        ram="15G"
    shell:
        """
        talon_create_GTF --db {input.db} \
            --whitelist {input.csv} \
            --build {params.ref_name} \
            -a {params.annotation_name} \
            --o {params.prefix} && \
        mv {params.prefix}_talon.gtf {output} 2> {log}
        """


###############################################################################
# intersect gtf with bam to keep only expressed genes and not those from the annotation
rule filter_gtf:
    input:
        gtf="results/{software}/{software}.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam"
    output:
        final=temp("results/{software}/{software}.{annot}.filtered.gtf")
    params:
        temp = "results/{software}/{software}.{annot}.EO.gtf"
    threads:1
    resources:
        ram="75G"
    log: "logs/only_seen_exons_{software}.{annot}.log"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        sed 's/*/./g' {input.gtf} | awk '$3 =="exon" && $7!="." && $5>$4' > {params.temp}
        bedtools intersect -u -s -split -a {params.temp} -b {input.bam} > {output.final} 2> {log}
        rm {params.temp}
        """

# Talon find some transcripts in 2 chr/other strand, remove the duplicate
rule correct_gtf:
    input:
        "results/{software}/{software}.{annot}.filtered.gtf"
    output:
        "results/{software}/{software}.{annot}.filtered_corrected.gtf"
    threads:1
    resources:
        ram="10G"
    script:
        "scripts/correct_gtf.py"

###############################################################################
# compare annotation to ref
rule gffcompare:
    input:
        test="results/{software}/{software}.{annot}.filtered_corrected.gtf",
        ref="results/utilities/uncompress.{annot}.gtf"
    output:
        result="results/gffcompare/{software}.{annot}.stats"
    threads:1
    resources:
        ram="15G"
    log: "logs/gffcompare_{software}.{annot}.log"
    conda:
        "envs/gffcompare.yaml"
    shadow: "shallow"
    params:
        prefix="{software}_{annot}",
    shell:
         """
         gffcompare {input.test} -r {input.ref} -o {params.prefix}
         mv {params.prefix}.stats {output.result} 2> {log}
         """


# format values to use it as R input
rule parse_gffcompare:
    input:
        expand("results/gffcompare/{software}.{annot}.stats",
                software=SOFTWARE, annot=config["annotation"])
    output:
        Sensitivity="results/gffcompare/Sensitivity.gffparse.tsv",
        Values="results/gffcompare/Values.gffparse.tsv"
    threads:1
    resources:
        ram="10G"
    log: "logs/parse_gffcompare.log"
    script:
        "scripts/gffcompare_parse.py"


# plot the expected results with ggplot2
rule gffcompare_report:
    input:
        Sensitivity="results/gffcompare/Sensitivity.gffparse.tsv",
        Values="results/gffcompare/Values.gffparse.tsv"
    output:
        "results/gffcompare/Graph.recap.pdf"
    threads:1
    resources:
        ram="10G"
    log: "logs/gffcompare_report.log"
    conda:
        "envs/r.yaml"
    resources:
        ram="6G"
    script:
        "scripts/gffcompare_report.R"


###########################################################################
rule install_SQANTI3:
    output:
        touch("results/utilities/SQANTI3_installed.txt")
    conda:
        "envs/sqanti.yaml"
    log: "logs/install_SQANTI3.log"
    threads:1
    resources:
        ram="5G"
    shell:
        """
        if ! [ -d "SQANTI3" ]; then
            echo "Download SQANTI3"
            git clone --depth 1 --branch v1.3 https://github.com/ConesaLab/SQANTI3.git
        fi
        if [ ! -f "./SQANTI3/utilities/gtfToGenePred" ]; then
            echo "Download gtfToGenePred"
            wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P ./SQANTI3/utilities/ -q
            echo "Add +x right to gtfToGenePred"
            chmod +x ./SQANTI3/utilities/gtfToGenePred
        fi
        if [ ! -d "./cDNA_Cupcake" ]; then
            echo "Download cDNA_Cupcake"
            git clone https://github.com/Magdoll/cDNA_Cupcake.git
            cd cDNA_Cupcake
            python setup.py build
            python setup.py install
            cd ..
        fi
        """


rule SQANTI3:
    input:
        isInstalled = "results/utilities/SQANTI3_installed.txt",
        gtfQuery = "results/{software}/{software}.{annot}.filtered_corrected.gtf",
        gtfRef = "results/utilities/uncompress.{annot}.gtf",
        fastaRef = 'results/utilities/reference.dna.uncompressed.fa',
    output:
        pdf = "results/SQANTI3/{software}/{software}.{annot}_sqanti_report.pdf",
        classification = "results/SQANTI3/{software}/{software}.{annot}_classification.txt",
        junction = "results/SQANTI3/{software}/{software}.{annot}_junctions.txt",
    params:
        dir="results/SQANTI3/{software}/",
        output_name="{software}.{annot}"
    conda:
        "envs/sqanti.yaml"
    log: "logs/SQANTI3_{software}.{annot}.log"
    threads:2
    resources:
        ram="16G"
    shell:
        """
        export PYTHONPATH=$PWD/cDNA_Cupcake/sequence/
        python ./SQANTI3/sqanti3_qc.py {input.gtfQuery} \
            {input.gtfRef} {input.fastaRef} \
            --gtf -d {params.dir} -o {params.output_name} 2> {log}
        """


rule parse_SQANTI3:
    input:
        classification = expand("results/SQANTI3/{software}/{software}.{annot}_classification.txt",
                software=SOFTWARE, annot=config["annotation"]),
        junctions = expand("results/SQANTI3/{software}/{software}.{annot}_junctions.txt",
                software=SOFTWARE, annot=config["annotation"])
    output:
        summary = "results/SQANTI3/summary.tsv"
    params:
        path = expand("results/SQANTI3/{software}/{software}.{annot}", software=SOFTWARE, annot=config["annotation"])
    conda:
        "envs/flair.yaml"
    log: "logs/parse_SQANTI3.log"
    threads:1
    resources:
        ram="8G"
    script:
        "scripts/sqanti_parse.py"


rule SQANTI_report:
    input:
        summary = "results/SQANTI3/summary.tsv"
    output:
        report = "results/SQANTI3/SQANTI_report.pdf"
    threads:1
    resources:
        ram="4G"
    log: "logs/SQANTI_report.log"
    conda:
        "envs/r.yaml"
    script:
        "scripts/sqanti_report.R"

################################################################################
rule RSEQC:
    input:
        bam = "results/minimap2/minimap.{annot}.sorted.bam",
        bed12 = "results/utilities/{annot}.converted.bed12"
    output:
        "results/RSeQC/{annot}.geneBodyCoverage.curves.pdf"
    threads:1
    resources:
        ram="8G"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        geneBody_coverage.py -i {input.bam} -r {input.bed12} -o results/RSeQC/{wildcards.annot}
        """


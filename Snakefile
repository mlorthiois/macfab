import glob

SOFTWARE=["stringtie", "flair", "talon", "bambu"]
wildcard_constraints:
    annot="[^./]+" # forbid wildcard "annot" to contain "/" or "." in order to ensure proper assignation

configfile: "config.yaml" # path to the config file
localrules: all # never launch the all and the compact rules on the cluster

 # a simple rule to launch the full pipeline without specifiying the final rule (graph)
rule all:
    input:
        results="results/gffcompare/Graph.recap.pdf",
        sqanti=expand("results/SQANTI3/{software}/{software}.{annot}_sqanti_report.pdf", 
                        annot=config["annotation"], 
                        software=SOFTWARE)


###############################################################################    
rule sam2bam:
    input:
        "results/minimap2/minimap.{annot}.sam"
    output:
        "results/minimap2/minimap.{annot}.sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads:10
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
        temp("results/utilities/reference.dna.uncompressed.fa")
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
    run:
        if input[0][-2:]=="gz":
            shell("gunzip -c {input} > {output}")
        else:
            shell("cp {input} {output}")


rule gtfToBed12:
    input:
        lambda wildcards: config["annotation"][wildcards.annot]
    output:
        "results/gtfToBed12/{annot}.converted.bed12"
    conda:
        "envs/minimap.yaml" # defines a conda env to launch the command
    threads:1
    shell:
        config['paftools_js'] + " gff2bed {input} > {output}"


###############################################################################
rule mapping:
    input:
        fastq=config["fastq_path"],
        bed="results/gtfToBed12/{annot}.converted.bed12",
        fa=config["reference_path"]
    output:
        "results/minimap2/minimap.{annot}.sam"
    conda:
        "envs/minimap.yaml"
    threads: workflow.cores
    shell:
        """
        minimap2 -t {threads} -ax splice --MD \
        --junc-bed {input.bed} {input.fa} {input.fastq} > {output}
        """


###############################################################################
# ensure that the env is set before bambu, avoid conflicts when two bambu instances are launched at the same time
rule install_bambu:
    output: touch("results/utilities/.R_config")
    conda:
        "envs/r.yaml"
    script:
        "scripts/install_bambu.R"


# asks for .R_config created by the R_config rule
rule bambu:
    input:
        gtf="results/utilities/uncompress.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam",
        fa="results/utilities/reference.dna.uncompressed.fa",
        isConfig="results/utilities/.R_config"
    output:
        bambu_gtf="results/bambu/bambu.{annot}.gtf"
    shadow: # snakemake symlinks the data, runs the soft inside the shadow dir then move the output out of it and delete it
        "shallow"
    conda:
        "envs/r.yaml"
    threads:1
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
    threads: workflow.cores
    shadow:
        "shallow"
    shell:
        """
        stringtie -L -G {input.gtf} -o {output.final_file} -p {threads} {input.bam}
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
        bedToGenePred {input} {params}
        genePredToGtf file {params} {output}
        """


###############################################################################
rule install_talon:
    output:
        touch("results/utilities/talon_installed.txt")
    conda:
        "envs/talon.yaml"
    shell:
        """
        if ! [ "$(command -v talon)" ]; then
            curl -L https://github.com/mortazavilab/TALON/archive/v5.0.zip -o ./results/utilities/TALON-5.0.zip
            cd results/utilities/
            unzip -q TALON-5.0.zip && rm TALON-5.0.zip
            cd TALON-5.0
            pip install .
            cd .. && rm -rf TALON-5.0
        fi
        """


rule talon_label_reads:
    input:
        bin="results/utilities/talon_installed.txt",
        sam="results/minimap2/minimap.{annot}.sam",
        fa=config["reference_path"]
    output:
        temp("results/talon/talon.{annot}_labeled.sam")
    params:
        prefix="results/talon/talon.{annot}"
    conda:
        "envs/talon.yaml"
    shadow:
        "shallow"
    shell:
        """
        talon_label_reads --f {input.sam} \
            --g {input.fa} --o {params.prefix}
        """


rule talon_initialize_database:
    input:
        bin="results/utilities/talon_installed.txt",
        gtf="results/utilities/uncompress.{annot}.gtf"
    output:
        db="results/talon/talon.{annot}.db"
    shadow:
        "shallow"
    params:
        reference_genome_name = config["reference_genome_name"],
        prefix = "results/talon/talon.{annot}",
        annotation_name = "{annot}"
    conda:
        "envs/talon.yaml"
    shell:
        """
        talon_initialize_database --f {input.gtf} \
            --g {params.reference_genome_name} --a {params.annotation_name} \
            --idprefix {params.prefix} --o {params.prefix}
        """


rule create_talon_configfile:
    input:
        "results/utilities/talon_installed.txt"
    params:
        id=config["sample_id"],
        description=config["sample_description"],
        ngs="Nanoseq",
        sam="results/talon/talon.{annot}_labeled.sam"
    shadow:
        "shallow"
    output:
        "results/talon/talon.{annot}.config"
    shell:
       "echo {params.id},{params.description},{params.ngs},{params.sam} > {output}" 


rule talon:
    input:
        bin="results/utilities/talon_installed.txt",
        config="results/talon/talon.{annot}.config",
        db="results/talon/talon.{annot}.db",
        sam="results/talon/talon.{annot}_labeled.sam"
    output:
        "results/talon/talon.{annot}_talon_read_annot.tsv"
    conda:
        "envs/talon.yaml"
    params:
        build=config["reference_genome_name"],
        prefix="results/talon/talon.{annot}"
    shadow:
        "shallow"
    shell:
        """
        talon --f {input.config} --db {input.db} \
            --build {params.build} --o {params.prefix}
        """


rule talon_create_GTF:
    input:
        bin="results/utilities/talon_installed.txt",
        db="results/talon/talon.{annot}.db"
    output:
        "results/talon/talon.{annot}.gtf"
    params:
        annotation_name="{annot}",
        ref_name=config["reference_genome_name"],
        prefix="results/talon/talon.{annot}"
    shadow:
        "shallow"
    conda:
        "envs/talon.yaml"
    shell:
        """
        talon_create_GTF --db {input.db} \
            -b {params.ref_name} -a {params.annotation_name} --o {params.prefix} \
            && mv {params.prefix}_talon.gtf {output}
        """


###############################################################################
# intersect gtf with bam to keep only expressed genes and not those from the annotation        
rule only_seen_exons:
    input:
        gtf="results/{software}/{software}.{annot}.gtf",
        bam="results/minimap2/minimap.{annot}.sorted.bam"
    output:
        final=temp("results/gffcompare/{software}.{annot}.filtered.gtf"),
        exon=temp("results/gffcompare/{software}.{annot}.EO.gtf")
    threads:1
    conda:
        "envs/bedtools.yaml"
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
    shell:
         """
         gffcompare {input.test} -r {input.ref} -o {params.prefix}
         mv {params.prefix}.stats {output.result}
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


###########################################################################
rule install_SQANTI3:
    output:
        touch("results/utilities/SQANTI3_installed.txt")
    conda:
        "envs/sqanti.yaml"
    shell:
        """
        if ! [ -d "SQANTI3" ]; then
            echo "Download SQANTI3"
            git clone https://github.com/ConesaLab/SQANTI3.git
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


# Remove lines where strand="." for SQANTI3
rule filter_strand_gtf:
    input:
        'results/{software}/{software}.{annot}.gtf'
    output:
        "results/{software}/{software}.{annot}_strand_corrected.gtf"
    shell:
        """
        awk '$7!="."' {input} > {output}
        """


rule SQANTI3:
    input:
        isInstalled = "results/utilities/SQANTI3_installed.txt",
        gtfQuery = "results/{software}/{software}.{annot}_strand_corrected.gtf",
        gtfRef = "results/utilities/uncompress.{annot}.gtf",
        fastaRef = 'results/utilities/reference.dna.uncompressed.fa',
    output:
        "results/SQANTI3/{software}/{software}.{annot}_sqanti_report.pdf"
    params:
        dir="results/SQANTI3/{software}/",
        output_name="{software}.{annot}"
    conda:
        "envs/sqanti.yaml"
    shell:
        """
        export PYTHONPATH=$PWD/cDNA_Cupcake/sequence/
        python ./SQANTI3/sqanti3_qc.py {input.gtfQuery} \
            {input.gtfRef} {input.fastaRef} \
            --gtf -d {params.dir} -o {params.output_name}
        """

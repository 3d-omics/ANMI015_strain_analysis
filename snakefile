######
# Metagenome variant callin pipeline
# Antton alberdi
# 2024/08/23
# Description: the pipeline creates creates a VCF file per sample.
#
# 1) Clone this repo.
# 2) Place mapping CRAMs in the input folder.
# 3) Create a screen session.
# 4) Launch the snakemake using the following code:
# module purge && module load snakemake/7.20.0 mamba/1.3.1
# snakemake -j 20 --cluster 'sbatch -o logs/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600
#
######

#List sample wildcards
SAMPLES, = glob_wildcards("input/{sample}.lib1.cram")

# CRAMs derived from 3D-omics/processed_data/MSEB0011-mg_quant/results/quantify/bowtie2/REF0009-mgg-pbdrep

# Define the reference genome
REFERENCE = "reference/REF0009-mgg-pbdrep.fa.gz"

# Define workflow
rule all:
    input:
        expand("output/{sample}.vcf.gz", sample=SAMPLES)

# Rule for converting CRAM to BAM
rule cram_to_bam:
    input:
        cram="input/{sample}.lib1.cram",
        ref=REFERENCE
    output:
        bam=temp("tmp/{sample}.bam")
    params:
        jobname="{sample}.1"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load samtools/1.15
        samtools view -b -T {input.ref} {input.cram} -o {output.bam}
        """

# Rule for sorting BAM
rule sort_bam:
    input:
        bam="tmp/{sample}.bam"
    output:
        sorted_bam=temp("tmp/{sample}.sorted.bam")
    params:
        jobname="{sample}.2"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load samtools/1.15
        samtools sort {input.bam} -o {output.sorted_bam}
        """

# Rule for indexing BAM
rule index_bam:
    input:
        sorted_bam="tmp/{sample}.sorted.bam"
    output:
        index=temp("tmp/{sample}.sorted.bam.bai")
    params:
        jobname="{sample}.3"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load samtools/1.15
        samtools index {input.sorted_bam}
        """

# Rule for adding and replacing read groups
rule fix_read_groups:
    input:
        sorted_bam="tmp/{sample}.sorted.bam",
        index="tmp/{sample}.sorted.bam.bai"
    output:
        fixed_bam=temp("tmp/{sample}.fixed.bam")
    params:
        jobname="{sample}.4"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load samtools/1.15
        samtools addreplacerg \
          -r 'ID:{wildcards.sample}_lib1\tLB:truseq_lib1\tPL:Illumina\tPU:unit1\tSM:{wildcards.sample}' \
          -o {output.fixed_bam} -w \
          {input.sorted_bam}
        """

# Rule for marking duplicates
rule mark_duplicates:
    input:
        fixed_bam="tmp/{sample}.fixed.bam"
    output:
        dedup_bam=temp("tmp/{sample}.dedup.bam"),
        metrics=temp("tmp/{sample}.dedup.txt")
    params:
        jobname="{sample}.5"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load jdk/1.8.0_291 picard/2.27.5
        picard MarkDuplicates I={input.fixed_bam} O={output.dedup_bam} M={output.metrics}
        """

# Rule for generating VCF
rule call_variants:
    input:
        fixed_bam="tmp/{sample}.dedup.bam",
        metrics="tmp/{sample}.dedup.txt",
        ref=REFERENCE
    output:
        vcf="output/{sample}.vcf.gz"
    params:
        jobname="{sample}.6"
    threads:
        1
    resources:
        mem_gb=4,
        time='00:05:00'
    shell:
        """
        module load bcftools/1.19
        bcftools mpileup -Ou -f {input.ref} {input.fixed_bam} | \
        bcftools call -mv -Oz -o {output.vcf} --ploidy 1
        rm {input.metrics}
        """

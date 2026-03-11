configfile: "configs/config.yaml"

EXPERIMENTS = config["samples"]


SAMPLES = {
    exp: glob_wildcards(
        f'{config["src_dir"]}/{exp}/fasta/{{sample}}_1.fastq'
    ).sample
    for exp in EXPERIMENTS
}



GENOMES = config["genomes"]



rule all:
    input:
        [
            expand(
                f'{config["src_dir"]}/{exp}/counts/{{genome}}_counts.csv',
                genome=GENOMES,
                sample=SAMPLES[exp]
            )
            for exp in EXPERIMENTS
        ]



rule STAR:
    input:
        r1 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_1.fastq',
        r2 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_2.fastq',
        index = f'{config["src_dir"]}/genomes/{{genome}}/index/'
    output:
        bam = f'{config["src_dir"]}/{{experiment}}/mapping/{{genome}}/{{sample}}ReadsPerGene.out.tab'
    params:
        tmpdir = lambda wildcards: f'{config["src_dir"]}/tmp/STAR_{wildcards.sample}_{wildcards.genome}'
    log:
        stderr="log/star/{experiment}/{genome}/{sample}.stderr"
    shell:
        """
        STAR \
            --runThreadN 8 \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --quantMode GeneCounts  \
            --genomeLoad NoSharedMemory \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir {params.tmpdir} \
            --outFileNamePrefix {config[src_dir]}/{wildcards.experiment}/mapping/{wildcards.genome}/{wildcards.sample}
            2> {log.stderr}
        """
rule preprocess_counts:
    input:
        folder = f'{config["src_dir"]}/{{experiment}}/mapping/{{genome}}'
    output:
        output = f'{config["src_dir"]}/{{experiment}}/counts/{{genome}}_counts.csv'

    shell:
        "Rscript /students/2025-2026/ros1_transcriptomics/minor_ROS1/scripts/R/preprocess_counts.R --input {input.folder} --output {output.output}"

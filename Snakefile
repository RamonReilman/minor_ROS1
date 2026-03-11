configfile: "configs/config.yaml"

EXPERIMENTS = config["samples"]
print(EXPERIMENTS)

SAMPLES = {
    exp: glob_wildcards(
        f'{config["src_dir"]}/{exp}/fasta/{{sample}}_1.fastq'
    ).sample
    for exp in EXPERIMENTS
}

print(SAMPLES)

GENOMES = config["genomes"]

print(GENOMES)

rule all:
    input:
        [
            expand(
                f'{config["src_dir"]}/{experiment}/mapping/{{genome}}/{{sample}}Aligned.sortedByCoord.out.bam',
                genome=GENOMES,
                sample=SAMPLES[experiment],
            )
            for experiment in EXPERIMENTS
        ]

rule STAR:
    input:
        r1 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_1.fastq',
        r2 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_2.fastq',
        index = f'{config["src_dir"]}/genomes/{{genome}}/index/'
    output:
        bam = f'{config["src_dir"]}/{{experiment}}/mapping/{{genome}}/{{sample}}Aligned.sortedByCoord.out.bam'
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

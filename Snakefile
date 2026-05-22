configfile: "configs/config_real_data.yaml"
try:
    EXPERIMENTS = config["samples"]
except:
    EXPERIMENTS = [""]

print(EXPERIMENTS)
SAMPLES = {
    exp: glob_wildcards(
        f'{config["src_dir"]}/{exp}/fasta/{{sample}}_R1_001.fastq'
    ).sample
    for exp in EXPERIMENTS
}



GENOMES = config["genomes"]


print(GENOMES, SAMPLES)
rule all:
    input:
        [
            expand(
                f'{config["src_dir"]}/{exp}/mapping/{{genome}}/{{sample}}ReadsPerGene.out.tab',
                genome=GENOMES,
                sample=SAMPLES[exp]
            )
            for exp in EXPERIMENTS
        ]



rule STAR:
    conda:
        "envs/STAR.yaml"
    input:
        r1 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_R1_001.fastq',
        r2 = f'{config["src_dir"]}/{{experiment}}/fasta/{{sample}}_R2_001.fastq',
        index = f'{config["src_dir"]}/genomes/{{genome}}/index/'
    output:
        bam = f'{config["src_dir"]}{{experiment}}/mapping/{{genome}}/{{sample}}ReadsPerGene.out.tab'
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
            --outFileNamePrefix {config[src_dir]}/{wildcards.experiment}/mapping/{wildcards.genome}/{wildcards.sample}
            2> {log.stderr}
        """

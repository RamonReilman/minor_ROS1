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

rule all:
    input:
        expand(
            "align/{experiment}/{genome}/{sample}.bam",
            experiment=EXPERIMENTS,
            genome=lambda wc: GENOMES[wc.experiment],
            sample=lambda wc: SAMPLES[wc.experiment]
        )

rule STAR:
    input:
        forward = "",
        reverse = ""
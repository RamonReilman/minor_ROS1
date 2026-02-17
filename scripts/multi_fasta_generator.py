import argparse
from pathlib import Path
from Bio import Entrez, SeqIO

def setup_parser():
    parser = argparse.ArgumentParser(
        prog = "Multi-Fasta generator",
        description = "Generate multi-fasta files based on gene name input",
        epilog = "help"
    )
    parser.add_argument("-i", "--input_list", required=True, help = "text file, with a gene names vertically seperated")
    parser.add_argument("-o", "--output", default="output.fasta", help = "Multifasta file output location")
    parser.add_argument("-e", "--email", default = "./config.ini", help = "Email used for Entrez", required=True)

    return parser

def validate_paths(paths):
    valid_paths = []
    for path in paths:
        valid_paths.append(Path(path))
        
    return valid_paths

def extract_gene_names(gene_file):
    gene_names = []
    with open(gene_file, "r") as f:
        for line in f:
            line = line.strip()
            gene_names.append(line)
    return gene_names

def get_fasta(gene_list):
    fasta = ""
    for gene in gene_list:
        handle = Entrez.esearch(db="nucleotide", term=f"{gene}[Gene Name] AND biomol_mrna[PROP] AND srcdb_refseq[PROP] AND human[ORGN]", retmax = 50)
        idlist = Entrez.read(handle)["IdList"]
        handle2 = Entrez.efetch(db="nucleotide", id=idlist, rettype="gb", retmode="text")
        try:
            header, seq = get_seqs(handle2)
        except Exception as _:
            print(f"{gene} not found, try manually")
            continue
        fasta += f">{gene} {header}\n{seq}\n"
    return fasta
        
def get_seqs(handle):
    seqs = list(SeqIO.parse(handle, "gb"))
    for seq in seqs:
        if "transcript variant 1," in seq.description.lower():
            return seq.id + "_" + seq.annotations["source"], seq.seq
    if "nm" in seqs[0].id.lower():
        return seqs[0].id + "_" + seqs[0].annotations["source"], seqs[0].seq
    raise ValueError()

def write_multifasta(fasta, output):
    with open(output, "w") as f:
        f.write(fasta)

def main():
    parser = setup_parser().parse_args()
    Entrez.email = parser.email
    
    input_file = parser.input_list
    output_file = parser.output
    input_file, output_file = validate_paths([input_file, output_file])
    gene_names = extract_gene_names(input_file)
    fasta = get_fasta(gene_names)
    write_multifasta(fasta, output_file)
    
    
    
    
    
    
    
    


if __name__ == "__main__":
    main()
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

input_fastq = snakemake.input[0]
output_csv = snakemake.output[0]  
output_fastq = snakemake.output[1]

data = []
trimmed_records = []

with open(input_fastq, 'r') as in_fq, open(output_fastq, 'w') as out_fq, open(output_csv, 'w') as out_csv:
    out_csv.write("ID,First_20_BP\n") # Write header to csv
    parser = FastqGeneralIterator(in_fq)
    for header, seq, qual in parser:
        first_20_bp = seq[:20]
        remaining_seq = seq[45:]
        remaining_quality = qual[45:]  # Correctly define remaining_quality

        # Write header and barcode to csv
        out_csv.write(f"{header},{first_20_bp}\n")

        # Write trimmed record to fastq
        out_fq.write(f"@{header}\n{remaining_seq}\n+\n{remaining_quality}\n")
import sys
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

input_fastq = snakemake.input[0]
output_csv = snakemake.output[0]  
output_fastq = snakemake.output[1]
upstream_len = snakemake.params[0]
bc_len = snakemake.params[1]
downstream_len = snakemake.params[2]

data = []
trimmed_records = []

with gzip.open(input_fastq, 'rt') as in_fq, gzip.open(output_fastq, 'wt') as out_fq, open(output_csv, 'w') as out_csv:
    out_csv.write("ID,Barcode\n") # Write header to csv
    parser = FastqGeneralIterator(in_fq)
    for header, seq, qual in parser:
        barcode = seq[upstream_len:upstream_len + bc_len]
        remaining_seq = seq[upstream_len + bc_len + downstream_len:]
        remaining_quality = qual[upstream_len + bc_len + downstream_len:]  # Correctly define remaining_quality

        # Write header and barcode to csv
        out_csv.write(f"{header},{barcode}\n")

        # Write trimmed record to fastq
        out_fq.write(f"@{header}\n{remaining_seq}\n+\n{remaining_quality}\n")

import pysam

def make_row(read, barcode):
    if any('MD' == tg[0] for tg in read.get_tags()):
        md = read.get_tag('MD')
        aligned = read.get_aligned_pairs(with_seq=True)
        indel_mm = []
        for x in aligned:
            if (None in x) or (x[2] in ['a', 'c', 't', 'g']):
                indel_mm.append(x)
    else:
        md = None
        indel_mm = None

    values = [
        read.query_name,
        read.seq,
        ''.join(chr(q + 33) for q in read.query_qualities),
        read.cigarstring,
        read.flag,
        read.reference_id,
        read.reference_start,
        read.reference_name,
        md,
        indel_mm,
        barcode + '\n'
    ]
    return ",".join(str(x) for x in values)

barcodes_file = snakemake.input[0]
alignment_file = snakemake.input[1]
merged_csv = snakemake.output[0]
samfile = pysam.AlignmentFile(alignment_file, 'r')

with open(barcodes_file, 'rt') as f_in, open(merged_csv, 'wt') as f_out:
    _ = next(f_in) # csv header
    f_out.write('Read ID,Sequence,Quality,CIGAR,Flag,Reference ID,Position,Probe_name,MD,Indel_MM,Barcode\n')

    for i, barcode_line in enumerate(f_in):
        bc_id, barcode = barcode_line.strip().split(',')
        bc_id = bc_id.split('/')[0]
        r1_aln = next(samfile)
        r2_aln = next(samfile)
        if bc_id != r1_aln.query_name or bc_id != r2_aln.query_name:
            raise ValueError(f"Headers don't match for {bc_id}!")

        f_out.write(make_row(r1_aln, barcode))
        f_out.write(make_row(r2_aln, barcode))

samfile.close()

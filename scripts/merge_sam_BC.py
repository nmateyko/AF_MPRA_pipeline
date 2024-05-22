#!/home/najmeh20/.conda/envs/GSC2024/bin/python
import pandas as pd
import pysam

def main():
    csv_file = '/scratch/st-cdeboer-1/najmeh/GSC-2024-04-12/pipeline/script/BWA/extracted_A7_R1_00_BC.csv'
    df_csv = pd.read_csv(csv_file)

    # Remove '/1' from the 'ID' column in the CSV
    df_csv['ID'] = df_csv['ID'].str.replace('/1$', '', regex=True)

    sam_file = '/scratch/st-cdeboer-1/najmeh/GSC-2024-04-12/pipeline/script/BWA/R1_R2_00.output.sam'
    samfile = pysam.AlignmentFile(sam_file, "r")

    data = []
    for read in samfile.fetch():
        data.append({
            'Read ID': read.query_name,
            'Sequence': read.seq,
            'Quality': ''.join(chr(q + 33) for q in read.query_qualities),
            'CIGAR': read.cigarstring,
            'Flag': read.flag,
            'Reference ID': read.reference_id,
            'Position': read.reference_start,
            'Probe_name': read.reference_name  #  reference name (third column in SAM)
        })

    df_sam = pd.DataFrame(data)

    samfile.close()

    # Merge df on 'Read ID'
    merged_df = pd.merge(df_csv, df_sam, left_on='ID', right_on='Read ID', how='inner')

    output_csv = '/scratch/st-cdeboer-1/najmeh/GSC-2024-04-12/pipeline/script/BWA/merged_R1_R2_00_sam_BC.csv'
    merged_df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    main()


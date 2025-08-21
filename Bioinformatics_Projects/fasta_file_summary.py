
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from io import StringIO

def parse_fasta(input_data):
    """
    Parse FASTA input, which can be a file path or a FASTA-formatted string.
    Handles multi-FASTA inputs.
    """
    if os.path.isfile(input_data):
        with open(input_data, 'r') as handle:
            return list(SeqIO.parse(handle, 'fasta'))
    else:
        return list(SeqIO.parse(StringIO(input_data), 'fasta'))

def compute_stats(record):
    """
    Compute statistics for a single SeqRecord.
    """
    return {
        'Sequence_ID': record.id,
        'Length': len(record.seq),
        'GC_Content': gc_fraction(record.seq) * 100,
        'Description': record.description
    }

def generate_report(input_data, output_csv=None, to_frame=True):
    """
    Generate a CSV report from FASTA input.
    """
    records = parse_fasta(input_data)
    data = [compute_stats(record) for record in records]
    df = pd.DataFrame(data)
    if to_frame:
        return df
    else:
        df.to_csv(output_csv, index=False)
        print(f"Report generated: {output_csv}")





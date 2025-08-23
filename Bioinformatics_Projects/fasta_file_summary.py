
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_fasta(*input_data):
    """
    Parse one or more FASTA inputs.

    Accepts:
    - file paths (str)
    - FASTA-formatted strings
    - SeqRecord objects
    - lists/tuples containing any of the above

    Returns a list of SeqRecord objects.
    """
    records = []

    # Normalize sources (flatten one level for lists/tuples)
    sources = []
    for item in input_data:
        if isinstance(item, (list, tuple)):
            sources.extend(item)
        else:
            sources.append(item)

    for src in sources:
        if isinstance(src, SeqRecord):
            records.append(src)
            continue

        if isinstance(src, str) and os.path.isfile(src):
            with open(src, 'r') as handle:
                records.extend(list(SeqIO.parse(handle, 'fasta')))
            continue

        if isinstance(src, str):
            stripped = src.strip()
            if stripped.startswith('>'):
                records.extend(list(SeqIO.parse(StringIO(src), 'fasta')))
            elif stripped:
                # Treat as raw sequence string without FASTA header
                seq_id = f"seq_{len(records) + 1}"
                records.append(SeqRecord(Seq(stripped), id=seq_id, description=""))
            continue

        raise TypeError(f"Unsupported input type: {type(src)}")

    return records

def compute_stats(*records):
    """
    Compute statistics for one or more SeqRecord objects.

    Returns:
    - dict if a single record is provided
    - list[dict] if multiple records are provided
    """
    # Flatten one level if a list/tuple was passed as a single arg
    if len(records) == 1 and isinstance(records[0], (list, tuple)):
        records = tuple(records[0])

    results = []
    for record in records:
        results.append({
            'Sequence_ID': record.id,
            'Length': len(record.seq),
            'GC_Content': gc_fraction(record.seq) * 100,
            'Description': record.description
        })

    if not results:
        return []
    return results[0] if len(results) == 1 else results

def generate_report(*input_data, output_csv=None, to_frame=True):
    """
    Generate a CSV report from one or more FASTA inputs.
    """
    records = parse_fasta(*input_data)
    stats = compute_stats(*records)
    data = [stats] if isinstance(stats, dict) else stats
    df = pd.DataFrame(data)
    if to_frame:
        return df
    else:
        if output_csv is None:
            raise ValueError("output_csv must be provided when to_frame is False.")
        df.to_csv(output_csv, index=False)
        print(f"Report generated: {output_csv}")





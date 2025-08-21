
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
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

def design_primer_pairs(input_data, 
                        primer_length=20, 
                        tm_min=55, 
                        tm_max=65, 
                        min_amplicon=100, 
                        max_amplicon=1000, 
                        top_n=10, 
                        max_candidates=50
                        ):
    """
    Design primer pairs for a given DNA sequence.
    Scans for candidate forward and reverse primers within Tm range,
    forms pairs within amplicon length range, and reports top N based on Tm similarity and average Tm.
    """
    records = parse_fasta(input_data)
    if not records:
        raise ValueError("No sequences found in input.")
    record = records[0]  # Take the first sequence
    seq = str(record.seq.upper())
    len_seq = len(seq)

    # Collect forward candidates
    forward_candidates = []
    for i in range(len_seq - primer_length + 1):
        primer_f = seq[i:i + primer_length]
        tm_f = MeltingTemp.Tm_NN(primer_f)
        if tm_min <= tm_f <= tm_max:
            forward_candidates.append({
                'start': i + 1,
                'seq': primer_f,
                'tm': tm_f
            })

    # Limit to top max_candidates by Tm
    if len(forward_candidates) > max_candidates:
        forward_candidates.sort(key=lambda x: -x['tm'])
        forward_candidates = forward_candidates[:max_candidates]

    # Collect reverse candidates
    rc = str(Seq(seq).reverse_complement())
    reverse_candidates = []
    for i in range(len_seq - primer_length + 1):
        primer_r = rc[i:i + primer_length]
        tm_r = MeltingTemp.Tm_NN(primer_r)
        if tm_min <= tm_r <= tm_max:
            start = len_seq - i - primer_length + 1
            reverse_candidates.append({
                'start': start,
                'seq': primer_r,
                'tm': tm_r
            })

    # Limit to top max_candidates by Tm
    if len(reverse_candidates) > max_candidates:
        reverse_candidates.sort(key=lambda x: -x['tm'])
        reverse_candidates = reverse_candidates[:max_candidates]

    # Generate valid pairs
    pairs = []
    for f in forward_candidates:
        for r in reverse_candidates:
            if r['start'] > f['start'] + primer_length - 1:
                amp_len = r['start'] + primer_length - f['start']
                if min_amplicon <= amp_len <= max_amplicon:
                    tm_diff = abs(f['tm'] - r['tm'])
                    avg_tm = (f['tm'] + r['tm']) / 2
                    pairs.append({
                        'forward_start': f['start'],
                        'forward_seq': f['seq'],
                        'forward_tm': f['tm'],
                        'reverse_start': r['start'],
                        'reverse_seq': r['seq'],
                        'reverse_tm': r['tm'],
                        'amplicon_length': amp_len,
                        'tm_diff': tm_diff,
                        'avg_tm': avg_tm
                    })

    # Sort pairs: smallest tm_diff first, then highest avg_tm
    pairs.sort(key=lambda p: (p['tm_diff'], -p['avg_tm']))

    # Take top N
    top_pairs = pairs[:top_n]

    return pd.DataFrame(top_pairs)

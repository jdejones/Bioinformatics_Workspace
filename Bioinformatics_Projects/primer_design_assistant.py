
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import primer3
from io import StringIO
from tqdm import tqdm

def _calc_tm_with_primer3(primer_sequence):
    """Return melting temperature (Tm) using primer3-py, with robust API handling."""
    seq_upper = str(primer_sequence).upper()
    # Try top-level function (preferred modern API)
    try:
        calc_tm_fn = getattr(primer3, 'calc_tm', None) or getattr(primer3, 'calcTm', None)
        if calc_tm_fn is not None:
            return float(calc_tm_fn(seq_upper))
    except Exception:
        pass

    # Try ThermoAnalysis class API (alternate API across versions)
    try:
        ThermoAnalysis = getattr(primer3, 'ThermoAnalysis', None)
        if ThermoAnalysis is None and hasattr(primer3, 'thermoanalysis'):
            ThermoAnalysis = getattr(primer3.thermoanalysis, 'ThermoAnalysis', None)
        if ThermoAnalysis is not None:
            analyzer = ThermoAnalysis()
            if hasattr(analyzer, 'calcTm'):
                return float(analyzer.calcTm(seq_upper))
            if hasattr(analyzer, 'calc_tm'):
                return float(analyzer.calc_tm(seq_upper))
    except Exception:
        pass

    raise RuntimeError("primer3-py Tm calculation is unavailable. Ensure primer3-py is installed and up to date.")

def parse_fasta(*input_data):
    """
    Parse one or more inputs into a flat list of SeqRecord objects.

    Accepts:
    - file paths (str)
    - FASTA-formatted strings
    - raw sequence strings (no header)
    - SeqRecord objects
    - lists/tuples of the above
    """
    records = []

    # Normalize and flatten a single level of lists/tuples
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
            if not stripped:
                continue
            if stripped.startswith('>'):
                records.extend(list(SeqIO.parse(StringIO(src), 'fasta')))
            else:
                # Treat as raw sequence string without FASTA header
                seq_id = f"seq_{len(records) + 1}"
                records.append(SeqRecord(Seq(stripped.upper()), id=seq_id, description=""))
            continue

        raise TypeError(f"Unsupported input type: {type(src)}")

    return records

def _design_primer_pairs_for_single_sequence(seq, seq_id,
                                            primer_length, 
                                            tm_min, 
                                            tm_max,
                                            min_amplicon, 
                                            max_amplicon,
                                            top_n, 
                                            max_candidates):
    """Worker: design primer pairs for a single sequence string."""
    seq = str(seq).upper()
    len_seq = len(seq)

    # Collect forward candidates
    forward_candidates = []
    for i in range(len_seq - primer_length + 1):
        primer_f = seq[i:i + primer_length]
        tm_f = _calc_tm_with_primer3(primer_f)
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
        tm_r = _calc_tm_with_primer3(primer_r)
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
                        'Sequence_ID': seq_id,
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
    return pairs[:top_n]


def design_primer_pairs(*input_data,
                        primer_length=20,
                        tm_min=55,
                        tm_max=65,
                        min_amplicon=100,
                        max_amplicon=1000,
                        top_n=10,
                        max_candidates=50
                        ):
    """
    Design primer pairs for one or more DNA sequences.
    Scans for candidate forward and reverse primers within Tm range,
    forms pairs within amplicon length range, and reports top N per sequence.
    When multiple sequences are provided, processing is parallelized across sequences.
    """
    records = parse_fasta(*input_data)
    if not records:
        raise ValueError("No sequences found in input.")

    # Single sequence: run inline
    if len(records) == 1:
        record = records[0]
        rows = _design_primer_pairs_for_single_sequence(
            str(record.seq), record.id,
            primer_length, tm_min, tm_max,
            min_amplicon, max_amplicon,
            top_n, max_candidates
        )
        return pd.DataFrame(rows)

    # Multiple sequences: parallel execution with progress bar
    start = time.perf_counter()
    rows_all = []
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                _design_primer_pairs_for_single_sequence,
                str(rec.seq), rec.id,
                primer_length, tm_min, tm_max,
                min_amplicon, max_amplicon,
                top_n, max_candidates
            )
            for rec in records
        ]

        with tqdm(total=len(futures), desc="Designing primers", unit="seq") as pbar:
            for fut in as_completed(futures):
                rows_all.extend(fut.result())
                pbar.update(1)

    elapsed = time.perf_counter() - start
    print(f"Primer design completed for {len(records)} sequences in {elapsed:.2f}s")

    return pd.DataFrame(rows_all)

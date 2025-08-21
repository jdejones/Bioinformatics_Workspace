import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def plot_gc_content_sliding_window(fasta_file, window_size=100, step=10):
    record = SeqIO.read(fasta_file, "fasta")
    seq = record.seq.upper()  # Ensure uppercase for GC calculation
    positions = []
    gc_values = []
    seq_len = len(seq)
    if seq_len < window_size:
        raise ValueError("Sequence length is shorter than window size.")
    for i in range(0, seq_len - window_size + 1, step):
        window = seq[i:i + window_size]
        gc = gc_fraction(window)
        pos = i + (window_size // 2)
        positions.append(pos)
        gc_values.append(gc)
    plt.figure(figsize=(10, 5))
    plt.plot(positions, gc_values)
    plt.xlabel("Position in sequence (bp)")
    plt.ylabel("GC content (%)")
    plt.title(f"Sliding Window GC Content\nWindow: {window_size} bp, Step: {step} bp")
    plt.grid(True)
    plt.show()

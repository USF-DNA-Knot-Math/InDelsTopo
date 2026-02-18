# utils.py
from collections import Counter
from Bio import SeqIO
import pandas as pd


def read_subsequence(chrom, start, end, direction='+', folder_path_fasta=None, dict_fast=None):
    """
    Read a subsequence from a plasmid FASTA file.

    Args:
        chrom (str): plasmid name
        start (int): start index (0-based)
        end (int): end index (exclusive)
        direction (str): '+' or '-' 
        folder_path_fasta (str): path to FASTA files
        dict_fast (dict): mapping plasmid name -> FASTA filename

    Returns:
        str: DNA sequence
    """
    if folder_path_fasta is None or dict_fast is None:
        raise ValueError("folder_path_fasta and dict_fast must be provided")

    fasta_path = folder_path_fasta + dict_fast[chrom]
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = record.seq
        if direction == '+':
            return str(seq[start:end])
        elif direction == '-':
            N = len(seq)
            return str(seq[N-end:N-start])
    raise ValueError("Invalid direction; must be '+' or '-'")

# Build table around R-loop start
# --------------------------
def produce_table(bed_df, folder_path_fasta, dict_fast, max_window_size=10,location='start'):
    """
    Produce a DataFrame with sequences around R-loop start sites.

    Args:
        bed_df (pd.DataFrame): BED dataframe
        folder_path_fasta (str): folder containing FASTA files
        dict_fast (dict): plasmid -> FASTA filename
        max_window_size (int): nucleotides upstream/downstream

    Returns:
        pd.DataFrame: table with 'sequence_pre' and 'sequence_post'
    """
    chrom = bed_df['chrom'][0]
    direction = bed_df['direction'][0]

    #Get plasmid length
    for record in SeqIO.parse(folder_path_fasta + dict_fast[chrom], "fasta"):
        N = len(record.seq)

    if direction == '+':
        Table = bed_df[['chrom', 'start', 'end']].copy()
    elif direction == '-':
        bed_df['new_start'] = N - bed_df['end']
        bed_df['new_end'] = N - bed_df['start']
        Table = bed_df[['chrom', 'new_start', 'new_end']].copy()
        Table.columns = ['chrom', 'start', 'end']
    else:
        raise ValueError("Invalid direction")
        
    if location=='start':
        # Windows around start
        Table['window_start'] = (Table.start - max_window_size).clip(lower=0)
        Table['window_end'] = (Table.start + max_window_size).clip(upper=N)
    elif location=='end':
        # Windows around end
        Table['window_start'] = (Table.end - max_window_size).clip(lower=0)
        Table['window_end'] = (Table.end + max_window_size).clip(upper=N)
    else:
        raise ValueError("Location needs to be start or end")
        
    # Read sequences
    Table['sequence_pre'] = Table.apply(
        lambda row: read_subsequence(row['chrom'], row['window_start'], row['start'],
                                     direction='+', folder_path_fasta=folder_path_fasta, dict_fast=dict_fast),
        axis=1)
    Table['sequence_post'] = Table.apply(
        lambda row: read_subsequence(row['chrom'], row['start'], row['window_end'],
                                     direction='+', folder_path_fasta=folder_path_fasta, dict_fast=dict_fast),
        axis=1)
    return Table

# K-mer processing
# --------------------------
def get_unique_words_freqs(Words):
    """
    Count unique words and their frequencies.

    Args:
        Words (list of str): list of DNA words

    Returns:
        tuple: (unique_words, frequencies)
    """
    counts = Counter(Words)
    return list(counts.keys()), list(counts.values())


def get_words_combined(Table, max_w=7):
    """
    Generate combined pre/post k-mers of varying lengths around R-loop start.

    Args:
        Table (pd.DataFrame): table with 'sequence_pre' and 'sequence_post'
        max_w (int): maximum k-mer length

    Returns:
        tuple: (unique_words, frequencies)
    """
    Words = []
    sizes = [(pre, post) for w in range(1, max_w+1)
                      for pre in range(w-1, w+1)
                      for post in range(w-1, w+1)]
    sizes = list(set(sizes))
    if (0,0) in sizes:
        sizes.remove((0,0))

    for pre, post in sizes:
        V_pre = Table['sequence_pre'].tolist()
        V_pre = [word[-pre:] if pre > 0 else '' for word in V_pre]

        V_post = Table['sequence_post'].tolist()
        V_post = [word[:post] if post > 0 else '' for word in V_post]

        V = [V_pre[i] + V_post[i] for i in range(len(V_pre))]
        Words += V

    return get_unique_words_freqs(Words)



"""
translatefastq.py - OPTIMIZED VERSION

Translates raw NGS .fastq files into amino acid frequency tables.

Inputs (via Snakemake config):
    mer              : number of amino acids per peptide
    startflank       : DNA sequence immediately upstream of the insert
    endflank         : DNA sequence immediately downstream of the insert
    PHD7             : bool -- True if using NEB PhD7 library (restricts codons)
    collapse         : bool -- True to collapse likely sequencing errors
    filter_min_reads : bool -- True to drop sequences with < MIN_READ_COUNT reads
                       before collapsing (optional, independent of collapse)

Output:
    CSV file with DNA sequences, AA translations, read counts, and
    normalized frequencies.

Originally created by Lindsey Brinton, University of Virginia, 2015.
Refactored and optimized for speed and correctness.
"""
import collections
import math
import statistics
import subprocess
import time
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Leader peptidase signal sequence (PSEX81 "GCTCAGCCGGCCATG" ) CCTTTAGTGGTACCTTTCTAT for NEB
PEPTIDASE_SEQ    = "CCTTTAGTGGTACCTTTCTAT"

# Endflank mutation detection (one-base variants upstream of endflank) for NEB "GGTGGAGG[ACG]"
ENDFLANK_MUT_PAT = "GGTGGAGG[ACG]"

# Amber stop codon (TAG) is read as glutamine in amber suppressor strains
# (TG1, ER2738 carry supE44). We map '_' -> 'Q' to reflect this biology.
AMBER_AA = "Q"

# Fixed error correction threshold (fraction of parent frequency)
FILTER_PERCENT = 0.02

# Minimum read count used by the optional pre-collapse filter
MIN_READ_COUNT = 5


CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "X", "TAG": "Z",  # Z = amber stop
    "TGC": "C", "TGT": "C", "TGA": "X", "TGG": "W",
    # Degenerate NNK codons
    "CCN": "P", "GGN": "G", "ACN": "T", "CGN": "R",
    "GCN": "A", "CTN": "L", "TCN": "S", "GTN": "V",
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def ascii_to_phred(char: str) -> int:
    """Convert an ASCII character to its Phred quality score."""
    return ord(char) - 33

def translate(seq: str) -> str:
    """
    Translate a DNA sequence to a protein sequence using CODON_TABLE.
    Codons not in the table are translated as 'X'.
    Amber stop codons ('Z') are replaced with AMBER_AA downstream.
    """
    if len(seq) % 3 != 0:
        return ""
    return "".join(
        CODON_TABLE.get(seq[i:i + 3], "X")
        for i in range(0, len(seq), 3)
    )

def grep_count(pattern: str, filepath: str, regex: bool = False) -> int:
    """Run grep -c (or -cE for regex) and return the integer count."""
    flag = "-cE" if regex else "-c"
    result = subprocess.run(
        ["grep", flag, pattern, filepath],
        capture_output=True, text=True, check=False
    )
    return int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0

def apply_min_reads_filter(
    freq_table: pd.DataFrame,
    min_count: int = MIN_READ_COUNT,
) -> pd.DataFrame:
    """
    Drop sequences with fewer than min_count reads.

    This is an optional pre-processing step that runs independently of
    collapsing. It removes noise sequences too sparse to be meaningful
    before any error-correction is applied.

    Parameters
    ----------
    freq_table : DataFrame sorted descending by count (column 0)
    min_count  : sequences with counts strictly below this value are dropped
    """
    original_count = len(freq_table)
    filtered = freq_table[freq_table.iloc[:, 0] >= min_count].copy()
    dropped = original_count - len(filtered)
    print(f"Min-reads filter: dropped {dropped} sequences with < {min_count} reads")
    print(f"Sequences remaining: {len(filtered)}")
    return filtered

def collapse_misreads(
    freq_table: pd.DataFrame,
    filter_percent: float = FILTER_PERCENT,
    top_n: int = None,
) -> pd.DataFrame:
    """
    Merge likely sequencing errors into their parent sequences.
    """
    # (top_n logic unchanged)
    if top_n is not None:
        freq_table = freq_table.head(top_n).copy()
        print(f"Keeping top {top_n} sequences: {len(freq_table)} retained")
    else:
        freq_table = freq_table.copy()

    seq_list  = list(freq_table.index)
    corrected = {seq: seq for seq in seq_list}

    counts = freq_table.iloc[:, 0].tolist()  # mutable, updated in-place

    thresholds = [
        filter_percent,
        filter_percent ** 2,
        filter_percent ** 3,
    ]

    total_collapsed = 0

    for idx in range(len(seq_list)):
        parent_seq  = seq_list[idx]
        parent_freq = counts[idx]   # reflects any counts merged into this parent
        parent_len  = len(parent_seq)

        if corrected[parent_seq] != parent_seq:
            continue

        if parent_freq * filter_percent < min(
            c for i, c in enumerate(counts)
            if corrected[seq_list[i]] == seq_list[i] and i > idx
        , default=0):
            print(
                f"Stopping at rank {idx + 1} "
                f"(parent count={parent_freq}, "
                f"2% threshold={parent_freq * filter_percent:.2f})"
            )
            break

        parent_collapsed = 0

        for i in range(idx + 1, len(seq_list)):
            seq = seq_list[i]

            if len(seq) != parent_len:
                continue
            if corrected[seq] != seq:
                continue

            diffs = sum(seq[j] != parent_seq[j] for j in range(parent_len))

            if 1 <= diffs <= 3 and counts[i] <= thresholds[diffs - 1] * parent_freq:
                corrected[seq]  = parent_seq
                counts[idx]    += counts[i]   # accumulate into parent immediately
                counts[i]       = 0            # zero out child so it can't be
                                               # mistaken for a parent later
                parent_freq     = counts[idx]  # keep local var in sync
                parent_collapsed += 1

        if parent_collapsed > 0:
            print(
                f"  Rank {idx + 1:>6} | count={parent_freq:>8} | "
                f"collapsed {parent_collapsed} sequence(s) into {parent_seq}"
            )
            total_collapsed += parent_collapsed

    print(f"\nTotal sequences collapsed: {total_collapsed}")
    print("Grouping corrected sequences...")
    freq_table.index = [corrected[seq] for seq in seq_list]
    freq_table = freq_table.groupby(freq_table.index).sum()

    return freq_table

def filter_non_nnk(freq_table: pd.DataFrame, mer: int) -> pd.DataFrame:
    """
    Remove reads containing codons not used by the NEB PhD7 NNK library.
    PhD7 uses NNK codons, so the third base of every codon must be G or T.
    """
    def is_nnk(seq: str) -> bool:
        return all(seq[i * 3 + 2] in ("G", "T") for i in range(mer))

    mask = [is_nnk(seq) for seq in freq_table.index]
    filtered = freq_table[mask]
    print(f"Sequences remaining after NNK filter: {len(filtered)}")
    return filtered

def export_to_csv(df: pd.DataFrame, filepath: str) -> None:
    """
    Export a DataFrame to CSV format (much faster than Excel, no row limits).

    CSV advantages:
    - 10-100x faster than Excel export
    - No 1M row limit per sheet
    - Smaller file size
    - Universal compatibility
    """
    if filepath.endswith('.xlsx'):
        filepath = filepath.replace('.xlsx', '.csv')

    df.to_csv(filepath, header=False)
    print(f"Exported to {filepath} ({len(df)} rows)")

# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def translatefastq(
    mer: int,
    filenamefastq: str,
    startflank: str,
    endflank: str,
    INSERTLESS_SEQ: str,
    filenameoutput: str,
    phd7: bool,
    collapse: bool,
    filter_min_reads: bool = False,
    top_n: int = None,
) -> pd.DataFrame:
    """
    Process a .fastq file from an NGS phage display screen.

    Steps
    -----
    1.  Count total reads, insertless controls, and peptidase signal reads.
    2.  Extract reads matching the flanking sequences.
    3.  Assess Phred quality scores.
    4.  Build frequency table.
    5.  (Optional) Drop sequences with < MIN_READ_COUNT reads.
    6.  (Optional) Collapse likely sequencing errors (1/2/3 bp, 2%/0.04%/0.0008%).
    7.  (Optional) Filter non-NNK codons for PhD7 libraries.
    8.  Translate DNA inserts to amino acid sequences.
    9.  Export results to CSV.

    Parameters
    ----------
    filter_min_reads : drop sequences with < MIN_READ_COUNT reads before any
                       collapsing. Runs even when collapse=False, so you can
                       use it as a standalone noise filter.
    top_n            : if set, retain only the top N sequences by frequency
                       before collapsing instead of processing all sequences.
                       Example: top_n=999999 keeps the top ~1M sequences.

    Returns
    -------
    pd.DataFrame with columns: count, AA, Normalized_Freq,
    5th_Percentile_Per_base_error
    """
    t_start = time.time()
    basepairs = 3 * mer
    print(f"\nImporting: {filenamefastq}")

    # ── Step 1: Read counts and QC metrics ───────────────────────────────────
    line_result  = subprocess.run(["wc", "-l", filenamefastq], capture_output=True, text=True)
    total_reads  = int(line_result.stdout.split()[0]) / 4
    print(f"Total reads: {total_reads:.0f}")

    insertless_count = grep_count(INSERTLESS_SEQ, filenamefastq)
    print(f"Insertless: {insertless_count / total_reads * 100:.2f}%")

    peptidase_count = grep_count(PEPTIDASE_SEQ, filenamefastq)
    print(f"With leader peptidase signal: {peptidase_count / total_reads * 100:.2f}%")

    misread_count = grep_count(f"[CGT]{PEPTIDASE_SEQ}", filenamefastq, regex=True)
    if peptidase_count > 0:
        print(f"Misread rate in peptidase signal: {misread_count / peptidase_count * 100:.2f}%")

    # ── Step 2: Extract reads with insert ────────────────────────────────────
    print("\nExtracting reads with inserts...")
    nt_pattern   = "[ACGT]"
    full_pattern = f"{startflank}({nt_pattern}){{{basepairs}}}{endflank}"
    print(f"Search pattern: {full_pattern}")

    grep_result = subprocess.run(
        ["grep", "-oE", full_pattern, filenamefastq],
        capture_output=True, text=True
    )
    if grep_result.returncode not in (0, 1):
        raise RuntimeError(f"grep failed: {grep_result.stderr}")

    raw_reads   = grep_result.stdout.splitlines()
    read_depth0 = len(raw_reads)
    pct_correct = read_depth0 / total_reads * 100
    print(f"Reads with insert: {read_depth0} ({pct_correct:.2f}%)")

    if pct_correct < 40:
        print(f"WARNING: Correct read percentage ({pct_correct:.2f}%) is below 40%.")

    # ── Step 3: NNK composition check ────────────────────────────────────────
    nnk_pattern = f"{startflank}([ACGT][ACGT][GT]){{{mer}}}{endflank}"
    nnk_count   = grep_count(nnk_pattern, filenamefastq, regex=True)
    print(f"NNK reads: {nnk_count / read_depth0 * 100:.2f}%")

    endflank_mut_count = grep_count(
        f"{startflank}.*{ENDFLANK_MUT_PAT}", filenamefastq, regex=True
    ) / 3
    mut_rate = endflank_mut_count / (endflank_mut_count + read_depth0)
    print(f"Endflank mutation rate: {mut_rate * 100:.4f}%")

    # ── Step 4: Phred quality scores ─────────────────────────────────────────
    print("\nCalculating Phred quality scores...")

    SAMPLE_SIZE = 10_000  # more than enough for a stable mean/std
    qual_scores_sampled = []

    with open(filenamefastq, "r") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq    = fh.readline().rstrip()
            plus   = fh.readline()
            qual   = fh.readline().rstrip()

            # Only use reads that contain the insert
            pos = seq.find(startflank)
            if pos == -1:
                continue

            # Slice the quality string to match the insert region
            q_insert = qual[pos : pos + len(startflank) + basepairs + len(endflank)]
            if len(q_insert) < len(startflank) + basepairs + len(endflank):
                continue

            qual_scores_sampled.append(q_insert)

            if len(qual_scores_sampled) >= SAMPLE_SIZE:
                break

    # Vectorized Phred calculation with numpy
    qual_array = np.array(
        [[ord(c) - 33 for c in qs] for qs in qual_scores_sampled],
        dtype=np.float32
    )
    per_read_mean     = qual_array.mean(axis=1)
    overall_avg_phred = per_read_mean.mean()
    phred_std         = per_read_mean.std()
    error_rate        = 10 ** (-overall_avg_phred / 10)
    adj_phred         = overall_avg_phred - 2 * phred_std
    adj_error_rate    = 10 ** (-adj_phred / 10)

    print(f"Sampled reads       : {len(qual_scores_sampled)}")
    print(f"Mean Phred score    : {overall_avg_phred:.2f}")
    print(f"Std dev             : {phred_std:.2f}")
    print(f"Mean error rate     : {error_rate:.6f}")
    print(f"5th pct error rate  : {adj_error_rate:.6f}")
    print(f"Filter threshold    : {FILTER_PERCENT * 100:.1f}%")

    if FILTER_PERCENT > adj_error_rate:
        print(
            f"WARNING: filter_percent ({FILTER_PERCENT}) exceeds 5th percentile "
            f"error rate ({adj_error_rate:.6f})."
        )

    # ── Step 5: Build frequency table ────────────────────────────────────────
    print("\nCounting reads...")
    counter    = collections.Counter(raw_reads)
    freq_table = (
        pd.DataFrame.from_dict(counter, orient="index")
        .sort_values(by=0, ascending=False)
    )
    print(f"Unique sequences: {len(freq_table)}")
    print(f"Total read depth: {freq_table.iloc[:, 0].sum()}")

    # ── Step 5.5: Optional min-reads filter ──────────────────────────────────
    if filter_min_reads:
        print(f"\nApplying min-reads filter (threshold: {MIN_READ_COUNT} reads)...")
        freq_table = apply_min_reads_filter(freq_table, MIN_READ_COUNT)

    # ── Step 6: Collapse misreads ─────────────────────────────────────────────
    if collapse:
        print("\nCollapsing misreads (1 bp: 2%, 2 bp: 0.04%, 3 bp: 0.0008%)...")
        freq_table = collapse_misreads(freq_table, FILTER_PERCENT, top_n)

    # ── Step 7: Filter non-NNK codons (PhD7 only) ────────────────────────────
    # Strip flanking sequences before codon-level filtering
    freq_table.index = freq_table.index.str[len(startflank):-len(endflank)]

    if phd7:
        print("\nFiltering non-NNK codons for PhD7 library...")
        freq_table = filter_non_nnk(freq_table, mer)

    # ── Step 8: Translate to amino acids ─────────────────────────────────────
    print("\nTranslating to amino acid sequences...")
    aa_sequences = [
        translate(seq).replace("_", AMBER_AA)
        for seq in freq_table.index
    ]

    freq_table["AA"]                            = aa_sequences
    freq_table["Normalized_Freq"]               = freq_table.iloc[:, 0] / read_depth0
    freq_table["5th_Percentile_Per_base_error"] = adj_error_rate

    result = freq_table.sort_values(by=0, ascending=False)
    print(f"\nFinal unique sequences : {len(result)}")
    print(f"Final total reads      : {result.iloc[:, 0].sum()}")

    # ── Step 9: Export ────────────────────────────────────────────────────────
    print("\nExporting results...")
    export_to_csv(result, filenameoutput)

    print(f"\nTotal time: {time.time() - t_start:.1f}s")
    return result

# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Inputs :", snakemake.input)
    print("Outputs:", snakemake.output)

    cfg = snakemake.config
    for i, (inp, out) in enumerate(zip(snakemake.input, snakemake.output)):
        print(f"\n--- Processing file {i + 1} of {len(snakemake.input)} ---")
        translatefastq(
            mer              = cfg["mer"],
            filenamefastq    = inp,
            startflank       = cfg["startflank"],
            endflank         = cfg["endflank"],
            INSERTLESS_SEQ   = cfg["insertless"],
            filenameoutput   = out,
            phd7             = cfg["PHD7"],
            collapse         = cfg["collapse"],
            filter_min_reads = cfg.get("filter_min_reads", False),
            top_n            = cfg.get("top_n", None),
        )

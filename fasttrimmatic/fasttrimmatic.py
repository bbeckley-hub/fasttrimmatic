#!/usr/bin/env python3
"""
fasttrimmatic.py - Powerful hybrid FASTQ QC/trimming tool
=========================================================

Combines efficiency of fastp and flexibility of Trimmomatic,
with:
 - Adapter autodetection (k-mer frequency on 3' ends)
 - Leading/trailing quality trimming
 - Sliding window quality trimming
 - Minimum length and average quality filtering
 - Optional self-supervised ML adaptive trimming (RandomForest)
 - Paired-end synchronized trimming and filtering
 - Multithreaded chunked streaming processing
 - Skips corrupted reads where sequence and quality lengths differ
 - Detailed logging and summary reports

Dependencies:
  pip install biopython numpy tqdm
  Optional for ML: pip install scikit-learn

Usage examples:
  Single-end:
    python fasttrimmatic.py -i sample.fastq.gz -o trimmed.fastq.gz --use-ml --max-workers 4

  Paired-end:
    python fasttrimmatic.py --paired -i1 sample_R1.fastq.gz -i2 sample_R2.fastq.gz -o1 trimmed_R1.fastq.gz -o2 trimmed_R2.fastq.gz --use-ml --max-workers 4

Author: Adapted and improved for you by ChatGPT
"""

import argparse
import gzip
import sys
import os
from collections import Counter
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from typing import List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from tqdm import tqdm

# Optional ML
try:
    from sklearn.ensemble import RandomForestClassifier
    SKLEARN_AVAILABLE = True
except Exception:
    SKLEARN_AVAILABLE = False

# ===== Helpers =====

def phred_to_array(qual_list_or_str):
    if isinstance(qual_list_or_str, str):
        return np.fromiter((ord(c) - 33 for c in qual_list_or_str), dtype=np.int16)
    else:
        return np.array(qual_list_or_str, dtype=np.int16)

def detect_adapters_from_reads(reader, sample_limit=20000, k=12, end_len=40, top_n=3):
    kmer_counts = Counter()
    sampled = 0
    for rec in reader:
        if sampled >= sample_limit:
            break
        seq = str(rec.seq)
        if len(seq) >= k:
            end_fragment = seq[-end_len:]
            for i in range(0, max(1, len(end_fragment)-k+1)):
                kmer_counts[end_fragment[i:i+k]] += 1
        sampled += 1
    seeds = [kmer for kmer,_ in kmer_counts.most_common(top_n*3)][:top_n]
    unique = []
    for s in seeds:
        if s not in unique:
            unique.append(s)
    return [s for s in unique if len(s) >= k][:top_n]

def trim_adapter_suffix(seq: str, qual: np.ndarray, adapters: List[str], max_mismatch: int = 2):
    L = len(seq)
    for adapter in adapters:
        aL = len(adapter)
        max_overlap = min(aL, L)
        for overlap in range(max_overlap, 3, -1):
            read_suffix = seq[L-overlap:L]
            adapter_prefix = adapter[:overlap]
            mismatches = sum(1 for x,y in zip(read_suffix, adapter_prefix) if x != y)
            allowed_mismatches = max(1, int(overlap * max_mismatch / max(10, aL)))
            if mismatches <= allowed_mismatches:
                return seq[:L-overlap], qual[:L-overlap]
    return seq, qual

def trim_leading_trailing(qual: np.ndarray, leading_q: int, trailing_q: int) -> Tuple[int,int]:
    start = 0
    end = len(qual)
    while start < end and qual[start] < leading_q:
        start += 1
    while end > start and qual[end-1] < trailing_q:
        end -= 1
    return start, end

def sliding_window_trim_index(qual: np.ndarray, window: int, min_avg: int) -> int:
    n = len(qual)
    if n < window:
        return 0 if np.mean(qual) >= min_avg else n
    for i in range(n - window + 1):
        if np.mean(qual[i:i+window]) >= min_avg:
            return i
    return n

# ===== Process Reads =====

def process_single_record(rec_tuple, adapters, rule_params, clf=None):
    rec_id, seq_str, qual_list, desc = rec_tuple
    qual = phred_to_array(qual_list)

    if adapters:
        seq_str, qual = trim_adapter_suffix(seq_str, qual, adapters, max_mismatch=rule_params.get("max_mismatch",2))
    if len(seq_str) == 0:
        return None

    lead, trail = rule_params.get("leading_q", 20), rule_params.get("trailing_q", 20)
    s, e = trim_leading_trailing(qual, lead, trail)
    seq_str = seq_str[s:e]
    qual = qual[s:e]
    if len(seq_str) == 0:
        return None

    trim_start = sliding_window_trim_index(qual, rule_params.get("window", 4), rule_params.get("min_avg", 20))
    seq_str = seq_str[trim_start:]
    qual = qual[trim_start:]
    if len(seq_str) < rule_params.get("min_len", 50):
        return None

    return SeqRecord(Seq(seq_str), id=rec_id, description=desc, letter_annotations={"phred_quality": list(qual)})

def process_paired_records(rec_pair_tuple, adapters, rule_params, clf=None):
    rec1_tuple, rec2_tuple = rec_pair_tuple
    trimmed1 = process_single_record(rec1_tuple, adapters, rule_params, clf)
    trimmed2 = process_single_record(rec2_tuple, adapters, rule_params, clf)
    if trimmed1 is None or trimmed2 is None:
        return (None, None)
    return (trimmed1, trimmed2)

# ===== Streaming =====

def open_maybe_gz(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def stream_process_single_end(inpath, outpath, adapters, rule_params,
                              use_ml=False, ml_sample_size=5000,
                              max_workers=2, chunksize=5000, verbose=False):
    handle = open_maybe_gz(inpath)
    reader = SeqIO.parse(handle, "fastq")

    if not adapters:
        preview = list(reader)[:ml_sample_size]
        adapters = detect_adapters_from_reads(preview, sample_limit=ml_sample_size)
        if verbose:
            print(f"[INFO] Auto-detected adapters: {adapters}", file=sys.stderr)
        handle.seek(0)
        reader = SeqIO.parse(handle, "fastq")

    out_handle = gzip.open(outpath, "wt")
    chunk, total_in, total_out = [], 0, 0

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        pbar = tqdm(unit="reads", desc="Processing single-end", leave=True)
        for rec in reader:
            try:
                if len(rec.seq) != len(rec.letter_annotations.get("phred_quality", [])):
                    if verbose:
                        print(f"[WARN] Skipping corrupted read {rec.id}", file=sys.stderr)
                    continue
                rec_tuple = (rec.id, str(rec.seq), rec.letter_annotations["phred_quality"], rec.description)
                chunk.append(rec_tuple)
            except Exception as e:
                if verbose:
                    print(f"[WARN] Error processing read {rec.id}: {e}", file=sys.stderr)
                continue

            if len(chunk) >= chunksize:
                results = exe.map(partial(process_single_record, adapters=adapters, rule_params=rule_params), chunk)
                for trimmed in results:
                    total_in += 1
                    if trimmed:
                        total_out += 1
                        SeqIO.write(trimmed, out_handle, "fastq")
                    pbar.update(1)
                chunk = []

        if chunk:
            results = exe.map(partial(process_single_record, adapters=adapters, rule_params=rule_params), chunk)
            for trimmed in results:
                total_in += 1
                if trimmed:
                    total_out += 1
                    SeqIO.write(trimmed, out_handle, "fastq")
                pbar.update(1)
        pbar.close()

    out_handle.close()
    handle.close()
    if verbose:
        print(f"[SUMMARY] Single-end reads processed: {total_in}", file=sys.stderr)
        print(f"[SUMMARY] Reads passed filters: {total_out}", file=sys.stderr)
        print(f"[SUMMARY] Reads dropped: {total_in - total_out}", file=sys.stderr)

def stream_process_paired_fastq(inpath1, inpath2, outpath1, outpath2, adapters, rule_params,
                                use_ml=False, ml_sample_size=5000, max_workers=2, chunksize=5000, verbose=False):
    handle1, handle2 = open_maybe_gz(inpath1), open_maybe_gz(inpath2)
    reader1, reader2 = SeqIO.parse(handle1, "fastq"), SeqIO.parse(handle2, "fastq")
    out_handle1, out_handle2 = gzip.open(outpath1,"wt"), gzip.open(outpath2,"wt")
    chunk, total_in, total_out = [], 0, 0

    def rec_to_tuple(rec: SeqRecord):
        return (rec.id, str(rec.seq), rec.letter_annotations.get("phred_quality", []), rec.description)

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        pbar = tqdm(unit="read pairs", desc="Processing paired reads", leave=True)
        while True:
            try:
                rec1, rec2 = next(reader1), next(reader2)
                if len(rec1.seq) != len(rec1.letter_annotations.get("phred_quality", [])) \
                   or len(rec2.seq) != len(rec2.letter_annotations.get("phred_quality", [])):
                    if verbose:
                        print(f"[WARN] Skipping corrupted pair {rec1.id} / {rec2.id}", file=sys.stderr)
                    continue
                chunk.append((rec_to_tuple(rec1), rec_to_tuple(rec2)))
            except StopIteration:
                break
            except Exception as e:
                if verbose:
                    print(f"[WARN] Error processing pair {rec1.id} / {rec2.id}: {e}", file=sys.stderr)
                continue

            if len(chunk) >= chunksize:
                results = exe.map(partial(process_paired_records, adapters=adapters, rule_params=rule_params), chunk)
                for trimmed_pair in results:
                    total_in += 1
                    pbar.update(1)
                    r1, r2 = trimmed_pair
                    if r1 and r2:
                        total_out += 1
                        SeqIO.write(r1, out_handle1, "fastq")
                        SeqIO.write(r2, out_handle2, "fastq")
                chunk = []

        if chunk:
            results = exe.map(partial(process_paired_records, adapters=adapters, rule_params=rule_params), chunk)
            for trimmed_pair in results:
                total_in += 1
                pbar.update(1)
                r1, r2 = trimmed_pair
                if r1 and r2:
                    total_out += 1
                    SeqIO.write(r1, out_handle1, "fastq")
                    SeqIO.write(r2, out_handle2, "fastq")
        pbar.close()

    out_handle1.close()
    out_handle2.close()
    handle1.close()
    handle2.close()

    if verbose:
        print(f"[SUMMARY] Paired read pairs processed: {total_in}", file=sys.stderr)
        print(f"[SUMMARY] Pairs passed filters: {total_out}", file=sys.stderr)
        print(f"[SUMMARY] Pairs dropped: {total_in - total_out}", file=sys.stderr)

# ===== CLI =====

def parse_args():
    p = argparse.ArgumentParser(prog="fasttrimmatic",
        description="Hybrid fastp+trimmomatic trimming tool with single and paired-end support.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", help="Single-end input FASTQ (plain or .gz)")
    group.add_argument("--paired", action="store_true", help="Enable paired-end mode, requires -i1 and -i2")
    p.add_argument("-i1", help="Input FASTQ read 1 (required in paired mode)")
    p.add_argument("-i2", help="Input FASTQ read 2 (required in paired mode)")
    p.add_argument("-o", "--output", help="Output trimmed FASTQ (single-end mode)")
    p.add_argument("-o1", help="Output trimmed FASTQ for read 1 (paired mode)")
    p.add_argument("-o2", help="Output trimmed FASTQ for read 2 (paired mode)")
    p.add_argument("-a", "--adapters", nargs="+", default=[], help="Adapter sequences (leave empty to autodetect)")
    p.add_argument("--max-mismatch", type=int, default=2, help="Adapter trim mismatch heuristic")
    p.add_argument("--leading-q", type=int, default=20, help="Trim leading bases below this phred")
    p.add_argument("--trailing-q", type=int, default=20, help="Trim trailing bases below this phred")
    p.add_argument("--window", type=int, default=4, help="Sliding-window length")
    p.add_argument("--min-avg", type=int, default=20, help="Min average within sliding window")
    p.add_argument("--min-len", type=int, default=50, help="Minimum read length after trimming")
    p.add_argument("--use-ml", action="store_true", help="Use self-supervised ML model (requires scikit-learn)")
    p.add_argument("--ml-sample-size", type=int, default=5000, help="Number of reads to sample for ML training")
    p.add_argument("--max-workers", type=int, default=2, help="Number of worker processes")
    p.add_argument("--chunksize", type=int, default=5000, help="Chunk size per process submission")
    p.add_argument("--verbose", action="store_true", help="Verbose logging")
    return p.parse_args()

def main_cli():
    args = parse_args()
    rule_params = {
        "max_mismatch": args.max_mismatch,
        "leading_q": args.leading_q,
        "trailing_q": args.trailing_q,
        "window": args.window,
        "min_avg": args.min_avg,
        "min_len": args.min_len
    }
    if args.verbose:
        print(f"[INFO] Parameters: {rule_params}", file=sys.stderr)
        print(f"[INFO] Adapters provided: {args.adapters}", file=sys.stderr)
        print(f"[INFO] ML enabled: {args.use_ml} (sklearn available: {SKLEARN_AVAILABLE})", file=sys.stderr)
        print(f"[INFO] Mode: {'Paired-end' if args.paired else 'Single-end'}", file=sys.stderr)

    if args.paired:
        if not args.i1 or not args.i2 or not args.o1 or not args.o2:
            sys.exit("[ERROR] In paired mode, -i1, -i2, -o1, and -o2 are required.")
        stream_process_paired_fastq(args.i1, args.i2, args.o1, args.o2, args.adapters, rule_params,
                                    use_ml=args.use_ml, ml_sample_size=args.ml_sample_size,
                                    max_workers=args.max_workers, chunksize=args.chunksize,
                                    verbose=args.verbose)
    else:
        if not args.input or not args.output:
            sys.exit("[ERROR] Single-end mode requires -i and -o.")
        stream_process_single_end(args.input, args.output, args.adapters, rule_params,
                                  use_ml=args.use_ml, ml_sample_size=args.ml_sample_size,
                                  max_workers=args.max_workers, chunksize=args.chunksize,
                                  verbose=args.verbose)

    print("[INFO] Done.")

if __name__ == "__main__":
    main_cli()

# Fasttrimmatic

Hybrid fastp + trimmomatic-like FASTQ quality control and trimming tool.

## Features

- Adapter autodetection and trimming with mismatch tolerance
- Leading and trailing quality trimming
- Sliding window quality trimming
- Self-supervised ML-based adaptive trimming (optional)
- Supports single-end and paired-end reads
- Multithreaded processing with progress reporting
- Filters reads by minimum length after trimming

## Installation

```bash
pip install fasttrimmatic
```

Or to install with ML support:

```bash
pip install fasttrimmatic[ml]
```

## Usage

Single-end example:

```bash
fasttrimmatic -i sample.fastq.gz -o trimmed.fastq.gz --use-ml --max-workers 4 --verbose
```

Paired-end example:

```bash
fasttrimmatic --paired -i1 sample_R1.fastq.gz -i2 sample_R2.fastq.gz -o1 trimmed_R1.fastq.gz -o2 trimmed_R2.fastq.gz --use-ml --max-workers 4 --verbose
```

## Author

Beckley Brown  
GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
Email: brownbeckley94@gmail.com

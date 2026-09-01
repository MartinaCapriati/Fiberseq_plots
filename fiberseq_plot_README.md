# Fiber-seq Condition Plotter

`fiberseq_plot_general.py` pools replicate Fiber-seq `sorted.nuc.bam` files by condition and creates:

- a stacked waterfall plot with one panel per condition
- a pooled metaplot with one line per condition
- one pooled matrix TSV per condition
- a small run summary for reproducibility

## Requirements

- Python 3.9+
- `pysmf` / PyFootPrint environment
- `numpy`
- `pandas`
- `matplotlib`
- `scipy` only when using `--cluster-waterfall hierarchical`



## Recommended Input Manifest

Create a tab-separated file such as `samples.tsv`:

```text
condition	sample	bam
untreated	rep1	/path/to/untreated_rep1.sorted.nuc.bam
untreated	rep2	/path/to/untreated_rep2.sorted.nuc.bam
treated	rep1	/path/to/treated_rep1.sorted.nuc.bam
treated	rep2	/path/to/treated_rep2.sorted.nuc.bam
```

## Example: Motif-Centered Plot

```bash
python fiberseq_plot_general.py \
  --ref /path/to/genome.fa \
  --motif chr11:98084116-98084126 \
  --flank 1000 \
  --manifest samples.tsv \
  --condition untreated \
  --condition treated \
  --outdir plots/example \
  --mark-motif BANP:chr11:98084116-98084126:red \
  --cluster-waterfall hierarchical \
  --cluster-window 100 \
  --smooth-window 5
```

## Example: Existing BED Region

```bash
python fiberseq_plot_general.py \
  --ref /path/to/genome.fa \
  --bed region.bed \
  --manifest samples.tsv \
  --outdir plots/example
```

## Inline BAM Input

For quick runs, skip the manifest and repeat `--bam`:

```bash
python fiberseq_plot_general.py \
  --ref /path/to/genome.fa \
  --bed region.bed \
  --bam untreated=rep1:/path/to/untreated_rep1.sorted.nuc.bam \
  --bam untreated=rep2:/path/to/untreated_rep2.sorted.nuc.bam \
  --bam treated=rep1:/path/to/treated_rep1.sorted.nuc.bam \
  --outdir plots/example
```

## Notes

- If `--condition` is omitted, condition order is inferred from the manifest or `--bam` arguments.
- If `--condition` is supplied, it controls plotting order and filters out unlisted conditions.
- For motif-centered runs, the script writes the generated BED file into the output directory.
- BAM indexes are not created automatically; the script warns if it cannot find one.

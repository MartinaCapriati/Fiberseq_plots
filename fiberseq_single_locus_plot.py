#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


LOGGER = logging.getLogger("fiberseq_plot")
plt = None
np = None
pd = None
linkage = None
leaves_list = None
pdist = None


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str = "region"

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end} ({self.name})"


@dataclass(frozen=True)
class Marker:
    name: str
    chrom: str
    start: int
    end: int
    color: str = "black"

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


@dataclass(frozen=True)
class SampleBam:
    condition: str
    sample: str
    bam: Path


def import_analysis_dependencies(require_scipy: bool = False) -> None:
    """Import heavy analysis dependencies only when the workflow runs."""
    global plt, np, pd, linkage, leaves_list, pdist

    if pd is None or np is None or plt is None:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as imported_plt
            import numpy as imported_np
            import pandas as imported_pd
        except ImportError as exc:
            raise ImportError(
                "Missing a required plotting dependency. Install matplotlib, numpy, "
                "and pandas in the environment used to run this script."
            ) from exc

        plt = imported_plt
        np = imported_np
        pd = imported_pd

    if require_scipy and (linkage is None or leaves_list is None or pdist is None):
        try:
            from scipy.cluster.hierarchy import leaves_list as imported_leaves_list
            from scipy.cluster.hierarchy import linkage as imported_linkage
            from scipy.spatial.distance import pdist as imported_pdist
        except ImportError as exc:
            raise ImportError("scipy is required for hierarchical waterfall clustering.") from exc

        linkage = imported_linkage
        leaves_list = imported_leaves_list
        pdist = imported_pdist


def parse_region_string(region_string: str, name: str = "region") -> Region:
    """Parse ``chrom:start-end`` strings, allowing comma separators."""
    region_string = region_string.strip()
    match = re.match(r"^([^:]+):([\d,]+)-([\d,]+)$", region_string)
    if not match:
        raise ValueError("Expected region format like chr7:55,443,600-55,444,100")

    start = int(match.group(2).replace(",", ""))
    end = int(match.group(3).replace(",", ""))
    if end <= start:
        raise ValueError(f"Region end must be greater than start: {region_string}")

    return Region(chrom=match.group(1), start=start, end=end, name=name)


def parse_marker_argument(marker_arg: str) -> Marker:
    """
    Parse a marker to show as a vertical line.

    Accepted formats:
        name:chrom:start-end:color
        name:chrom:start-end
        chrom:start-end:color
        chrom:start-end
    """
    parts = marker_arg.split(":")

    if len(parts) == 2:
        name = None
        chrom, coords = parts
        color = "black"
    elif len(parts) == 3:
        if "-" in parts[1]:
            name = None
            chrom, coords, color = parts
        else:
            name, chrom, coords = parts
            color = "black"
    elif len(parts) == 4:
        name, chrom, coords, color = parts
    else:
        raise ValueError(
            f"Invalid marker format: {marker_arg}. Use name:chrom:start-end:color, "
            "name:chrom:start-end, chrom:start-end:color, or chrom:start-end."
        )

    try:
        start_s, end_s = coords.split("-", 1)
        start = int(start_s.replace(",", ""))
        end = int(end_s.replace(",", ""))
    except ValueError as exc:
        raise ValueError(f"Invalid marker coordinates: {marker_arg}") from exc

    if end <= start:
        raise ValueError(f"Marker end must be greater than start: {marker_arg}")

    marker_name = name if name is not None else f"{chrom}:{start}-{end}"
    return Marker(name=marker_name, chrom=chrom, start=start, end=end, color=color)


def parse_marker_arguments(marker_args: Iterable[str] | None) -> list[Marker]:
    if not marker_args:
        return []
    return [parse_marker_argument(marker_arg) for marker_arg in marker_args]


def draw_motif_markers(ax: plt.Axes, markers: Sequence[Marker], y_text: bool = False) -> None:
    """Draw marker center lines on an axis."""
    for marker in markers:
        ax.axvline(
            marker.center,
            linestyle="--",
            linewidth=1.2,
            color=marker.color,
            alpha=0.85,
        )
        if y_text:
            ylim = ax.get_ylim()
            ax.text(
                marker.center,
                ylim[1],
                marker.name,
                color=marker.color,
                rotation=90,
                ha="right",
                va="top",
                fontsize=8,
            )


def write_motif_centered_bed(motif: str, flank: int, out_bed: Path) -> tuple[Path, Region, Region]:
    """Create a single-region BED centered on the supplied motif interval."""
    if flank < 0:
        raise ValueError("--flank must be >= 0")

    motif_region = parse_region_string(motif, name="motif")
    bed_region = Region(
        chrom=motif_region.chrom,
        start=max(0, motif_region.center - flank),
        end=motif_region.center + flank,
        name="motif_centered",
    )

    with out_bed.open("w", encoding="utf-8") as handle:
        handle.write(f"{bed_region.chrom}\t{bed_region.start}\t{bed_region.end}\t{bed_region.name}\n")

    return out_bed, motif_region, bed_region


def read_first_bed_region(bed: Path) -> Region:
    """Read the first non-comment BED interval."""
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"Invalid BED line: {line}")
            name = fields[3] if len(fields) > 3 else "region"
            return Region(fields[0], int(fields[1]), int(fields[2]), name)

    raise ValueError(f"No valid BED region found in {bed}")


def parse_condition_bam_argument(arg: str) -> SampleBam:
    """
    Parse ``--bam CONDITION=SAMPLE:/path.bam`` or ``--bam CONDITION=/path.bam``.
    """
    if "=" not in arg:
        raise ValueError(f"Invalid --bam: {arg}. Use CONDITION=SAMPLE:/path.bam")

    condition, rest = arg.split("=", 1)
    condition = condition.strip()
    if not condition:
        raise ValueError(f"Missing condition in --bam: {arg}")

    if ":" in rest:
        sample, bam_s = rest.split(":", 1)
        sample = sample.strip()
    else:
        bam_s = rest
        sample = derive_sample_name(Path(bam_s))

    bam = Path(bam_s).expanduser()
    if not sample:
        raise ValueError(f"Missing sample name in --bam: {arg}")
    if not bam.exists():
        raise FileNotFoundError(f"BAM not found: {bam}")

    return SampleBam(condition=condition, sample=sample, bam=bam)


def derive_sample_name(bam: Path) -> str:
    """Create a readable sample name from a BAM filename."""
    name = bam.name
    for suffix in (".sorted.nuc.bam", ".nuc.bam", ".bam"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return bam.stem


def parse_manifest_sep(sep: str) -> str:
    """Allow friendly escaped separators such as '\\t'."""
    if sep == r"\t":
        return "\t"
    return sep


def load_manifest(path: Path, sep: str) -> list[SampleBam]:
    """Read a sample manifest with condition, sample, and bam columns."""
    rows: list[SampleBam] = []
    delimiter = parse_manifest_sep(sep)

    with path.expanduser().open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")

        required = {"condition", "sample", "bam"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Manifest is missing required column(s): {', '.join(sorted(missing))}. "
                "Expected columns: condition, sample, bam."
            )

        for line_number, row in enumerate(reader, start=2):
            condition = (row.get("condition") or "").strip()
            sample = (row.get("sample") or "").strip()
            bam_text = (row.get("bam") or "").strip()
            if not condition or not sample or not bam_text:
                raise ValueError(f"Blank condition, sample, or bam at manifest line {line_number}")
            bam = Path(bam_text).expanduser()
            if not bam.exists():
                raise FileNotFoundError(f"BAM not found at manifest line {line_number}: {bam}")
            rows.append(SampleBam(condition=condition, sample=sample, bam=bam))

    if not rows:
        raise ValueError(f"Manifest contains no samples: {path}")

    return rows


def group_samples_by_condition(
    samples: Sequence[SampleBam],
    requested_conditions: Sequence[str] | None,
) -> "OrderedDict[str, list[SampleBam]]":
    """Group samples and either infer or validate condition order."""
    if not samples:
        raise ValueError("No BAMs supplied. Use --bam or --manifest.")

    seen_conditions: list[str] = []
    for sample in samples:
        if sample.condition not in seen_conditions:
            seen_conditions.append(sample.condition)

    condition_order = list(requested_conditions or seen_conditions)
    unknown = [condition for condition in condition_order if condition not in seen_conditions]
    if unknown:
        raise ValueError(
            f"Condition(s) requested but not found in BAM inputs: {', '.join(unknown)}. "
            f"Available: {', '.join(seen_conditions)}"
        )

    grouped: "OrderedDict[str, list[SampleBam]]" = OrderedDict((condition, []) for condition in condition_order)
    ignored = sorted(set(seen_conditions).difference(condition_order))
    if ignored:
        LOGGER.warning("Ignoring condition(s) not listed with --condition: %s", ", ".join(ignored))

    for sample in samples:
        if sample.condition in grouped:
            grouped[sample.condition].append(sample)

    empty = [condition for condition, condition_samples in grouped.items() if not condition_samples]
    if empty:
        raise ValueError(f"No BAMs supplied for condition(s): {', '.join(empty)}")

    return grouped


def ensure_bam_index(bam: Path) -> None:
    """Warn if the expected BAM index is missing."""
    bam_s = str(bam)
    bai1 = Path(bam_s + ".bai")
    bai2 = Path(bam_s[:-4] + ".bai") if bam_s.endswith(".bam") else Path(bam_s + ".bai")
    if not bai1.exists() and not bai2.exists():
        LOGGER.warning("No BAM index found for %s. Try: samtools index %s", bam, bam)


def import_pysmf():
    """Import pysmf only when matrix extraction is actually requested."""
    try:
        import pysmf  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Could not import pysmf. Activate the PyFootPrint/pysmf environment "
            "before running matrix extraction."
        ) from exc
    return pysmf


def load_project_compatible(bam: Path, ref: Path):
    """Call pysmf.load_project across common pysmf signatures."""
    pysmf = import_pysmf()
    bam_s = str(bam)
    ref_s = str(ref)

    try:
        return pysmf.load_project(bam_s, ref_s, "Fiberseq")
    except TypeError:
        try:
            return pysmf.load_project(bam_s, reference=ref_s, sequencing_type="Fiberseq")
        except TypeError:
            return pysmf.load_project(bam_s, reference_genome=ref_s, sequencing_type="Fiberseq")


def extract_matrix(
    bam: Path,
    ref: Path,
    bed: Path,
    region_index: int = 0,
    padding: int = 0,
    min_reads: int = 20,
) -> pd.DataFrame:
    """Extract a methylation matrix from one BAM using pysmf."""
    ensure_bam_index(bam)
    smf = load_project_compatible(bam, ref)
    smf.load_regions(str(bed))
    matrix = smf.extract_methylation_matrix(region_index, padding, min_reads)

    if matrix is None:
        raise ValueError(f"pysmf returned None for {bam}")
    if hasattr(matrix, "data_frame"):
        df = matrix.data_frame.copy()
    elif hasattr(matrix, "dataframe"):
        df = matrix.dataframe.copy()
    elif isinstance(matrix, pd.DataFrame):
        df = matrix.copy()
    else:
        raise TypeError(f"Unsupported extract_methylation_matrix return type: {type(matrix)}")

    cols = []
    for col in df.columns:
        try:
            cols.append(int(col))
        except (TypeError, ValueError):
            cols.append(col)
    df.columns = cols
    return sort_matrix_columns(df)


def sort_matrix_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Sort genomic-position columns while leaving mixed-column matrices unchanged."""
    numeric = []
    other = []
    for col in df.columns:
        try:
            numeric.append(int(col))
        except (TypeError, ValueError):
            other.append(col)

    if numeric and not other:
        return df.loc[:, sorted(numeric)]
    return df


def matrix_average_profile(df: pd.DataFrame) -> pd.Series:
    """Return the pooled average profile used by the metaplot."""
    df = sort_matrix_columns(df)
    numeric = df.apply(pd.to_numeric, errors="coerce")
    return 1 - numeric.mean(axis=0)


def smooth_series(series: pd.Series, window: int | None) -> pd.Series:
    """Centered rolling mean smoothing."""
    if window is None or window <= 1:
        return series
    return series.rolling(window=window, center=True, min_periods=1).mean()


def cluster_reads_by_center_accessibility(
    df: pd.DataFrame,
    center: int,
    window: int,
    method: str,
) -> pd.DataFrame:
    """Order waterfall reads around a center window."""
    if method == "none":
        return df

    df = sort_matrix_columns(df)
    numeric = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if numeric.shape[0] < 2:
        return df

    positions = np.array(numeric.columns, dtype=float)
    keep = (positions >= center - window) & (positions <= center + window)
    if keep.sum() < 2:
        LOGGER.warning("Too few columns within center +/- %s; skipping clustering.", window)
        return df

    cluster_df = numeric.loc[:, keep]

    if method == "score":
        scores = cluster_df.mean(axis=1)
        return df.loc[scores.sort_values(ascending=False).index]

    if method == "hierarchical":
        if linkage is None or leaves_list is None or pdist is None:
            raise ImportError("scipy is required for hierarchical waterfall clustering.")
        distances = pdist(cluster_df.values, metric="euclidean")
        if np.all(distances == 0):
            return df
        order = leaves_list(linkage(distances, method="average"))
        return df.iloc[order]

    raise ValueError(f"Unknown clustering method: {method}")


def plot_stacked_waterfalls(
    condition_matrices: Mapping[str, pd.DataFrame],
    out_path: Path,
    title: str,
    center: int,
    cluster_window: int,
    cluster_method: str,
    max_reads: int | None = None,
    motif_markers: Sequence[Marker] | None = None,
    random_seed: int = 1,
    dpi: int = 200,
) -> None:
    """Write a stacked waterfall figure."""
    motif_markers = motif_markers or []
    processed: "OrderedDict[str, pd.DataFrame]" = OrderedDict()
    heights: list[int] = []

    for condition, df in condition_matrices.items():
        numeric = sort_matrix_columns(df).apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")

        if max_reads is not None and numeric.shape[0] > max_reads:
            numeric = numeric.sample(n=max_reads, random_state=random_seed)

        numeric = cluster_reads_by_center_accessibility(
            numeric,
            center=center,
            window=cluster_window,
            method=cluster_method,
        )

        processed[condition] = numeric
        heights.append(max(1, numeric.shape[0]))

    n_conditions = len(processed)
    fig_height = max(6, min(32, sum(heights) / 45 + n_conditions * 1.5))
    fig, axes = plt.subplots(
        n_conditions,
        1,
        figsize=(14, fig_height),
        sharex=True,
        gridspec_kw={"height_ratios": heights},
    )

    if n_conditions == 1:
        axes = [axes]

    for ax, (condition, df) in zip(axes, processed.items()):
        if df.empty:
            ax.text(0.5, 0.5, f"{condition}: empty", ha="center", va="center")
            ax.set_ylabel(f"{condition}\nreads")
            continue

        x = np.array(df.columns, dtype=float)
        for y_idx, (_, row) in enumerate(df.iterrows()):
            values = row.to_numpy(dtype=float)
            present = values > 0
            absent = values <= 0

            if absent.any():
                ax.scatter(
                    x[absent],
                    np.full(absent.sum(), y_idx),
                    marker="|",
                    s=7,
                    alpha=0.30,
                    color="lightgrey",
                    linewidths=0.5,
                )
            if present.any():
                ax.scatter(
                    x[present],
                    np.full(present.sum(), y_idx),
                    marker="s",
                    s=4,
                    alpha=0.85,
                    color="tab:blue",
                    linewidths=0,
                )

        ax.axvline(center, linestyle="--", linewidth=1, color="black", alpha=0.6)
        draw_motif_markers(ax, motif_markers, y_text=False)
        ax.set_ylabel(f"{condition}\nreads")
        ax.set_ylim(-1, df.shape[0] + 1)

    axes[-1].set_xlabel("Genomic position")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_condition_metaplot(
    condition_matrices: Mapping[str, pd.DataFrame],
    out_path: Path,
    title: str,
    smooth_window: int,
    motif_markers: Sequence[Marker] | None = None,
    ylabel: str = "1 - mean matrix value",
    dpi: int = 200,
) -> None:
    """Write one pooled average profile per condition."""
    motif_markers = motif_markers or []

    fig, ax = plt.subplots(figsize=(12, 4.8))
    for condition, df in condition_matrices.items():
        profile = smooth_series(matrix_average_profile(df), smooth_window)
        ax.plot(
            np.array(profile.index, dtype=float),
            profile.to_numpy(dtype=float),
            linewidth=1.8,
            label=condition,
        )

    draw_motif_markers(ax, motif_markers, y_text=True)
    ax.set_title(title)
    ax.set_xlabel("Genomic position")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def collect_samples(args: argparse.Namespace) -> list[SampleBam]:
    """Load BAM definitions from CLI arguments and/or a manifest."""
    samples: list[SampleBam] = []

    if args.manifest:
        samples.extend(load_manifest(Path(args.manifest), args.manifest_sep))

    for bam_arg in args.bam or []:
        samples.append(parse_condition_bam_argument(bam_arg))

    return samples


def write_run_summary(
    out_path: Path,
    args: argparse.Namespace,
    region: Region,
    cluster_center: int,
    samples_by_condition: Mapping[str, Sequence[SampleBam]],
) -> None:
    """Write a small, human-readable summary for reproducibility."""
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("Fiber-seq plotting run\n")
        handle.write("======================\n\n")
        handle.write(f"Reference: {Path(args.ref).expanduser()}\n")
        handle.write(f"Region: {region.label}\n")
        handle.write(f"Output directory: {Path(args.outdir).expanduser()}\n")
        handle.write(f"Waterfall clustering: {args.cluster_waterfall}\n")
        handle.write(f"Cluster center: {cluster_center}\n")
        handle.write(f"Cluster window: {args.cluster_window}\n")
        handle.write(f"Minimum reads: {args.min_reads}\n")
        handle.write(f"Smoothing window: {args.smooth_window}\n\n")
        handle.write("Samples\n")
        handle.write("-------\n")
        for condition, samples in samples_by_condition.items():
            for sample in samples:
                handle.write(f"{condition}\t{sample.sample}\t{sample.bam}\n")


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be > 0")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Pool Fiber-seq replicates per condition and plot waterfalls/metaplots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--ref", required=True, help="Reference FASTA used by pysmf.")
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--bed", help="BED file containing the region to plot. The first BED record is used.")
    region_group.add_argument("--motif", help="Motif interval to center on, e.g. chr7:55443890-55443899.")
    parser.add_argument("--flank", type=nonnegative_int, default=5000, help="Flank in bp around --motif center.")

    parser.add_argument(
        "--manifest",
        help="Sample manifest with columns: condition, sample, bam. Can be combined with repeated --bam.",
    )
    parser.add_argument("--manifest-sep", default=r"\t", help="Manifest delimiter, e.g. '\\t' or ','.")
    parser.add_argument(
        "--bam",
        action="append",
        default=[],
        help="BAM definition. Use CONDITION=SAMPLE:/path.bam or CONDITION=/path.bam. Repeatable.",
    )
    parser.add_argument(
        "--condition",
        action="append",
        default=[],
        help="Condition order to plot. If omitted, order is inferred from --manifest/--bam input.",
    )

    parser.add_argument("--outdir", default="fiberseq_plots", help="Output directory.")
    parser.add_argument("--output-prefix", default="pooled", help="Prefix for figure filenames.")
    parser.add_argument("--region-index", type=nonnegative_int, default=0, help="Region index for pysmf extraction.")
    parser.add_argument("--padding", type=nonnegative_int, default=0, help="Padding passed to pysmf extraction.")
    parser.add_argument("--min-reads", type=nonnegative_int, default=20, help="Minimum reads passed to pysmf.")
    parser.add_argument("--smooth-window", type=positive_int, default=1, help="Rolling mean window for metaplot.")

    parser.add_argument(
        "--cluster-waterfall",
        choices=["none", "hierarchical", "score"],
        default="none",
        help="Read ordering for waterfall panels.",
    )
    parser.add_argument("--cluster-window", type=nonnegative_int, default=1000, help="Window around cluster center in bp.")
    parser.add_argument("--cluster-center", type=int, default=None, help="Genomic coordinate for clustering.")
    parser.add_argument("--max-waterfall-reads", type=positive_int, default=None, help="Subsample reads per panel.")
    parser.add_argument("--random-seed", type=int, default=1, help="Seed used when subsampling waterfall reads.")

    parser.add_argument(
        "--mark-motif",
        action="append",
        default=[],
        help=(
            "Vertical marker to draw. Repeatable. Formats: name:chrom:start-end:color, "
            "name:chrom:start-end, chrom:start-end:color, or chrom:start-end."
        ),
    )
    parser.add_argument("--waterfall-format", choices=["png", "pdf", "svg"], default="png")
    parser.add_argument("--metaplot-format", choices=["pdf", "png", "svg"], default="pdf")
    parser.add_argument("--dpi", type=positive_int, default=200, help="Figure DPI for raster outputs.")
    parser.add_argument("--metaplot-ylabel", default="1 - mean matrix value", help="Y-axis label for the metaplot.")
    parser.add_argument("--verbose", action="store_true", help="Print additional progress details.")

    return parser


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(format="%(levelname)s: %(message)s", level=level)


def run(args: argparse.Namespace) -> dict[str, Path]:
    """Run the workflow and return key output paths."""
    ref = Path(args.ref).expanduser()
    if not ref.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref}")

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    motif_center = None
    if args.motif:
        bed, motif_region, plot_region = write_motif_centered_bed(
            args.motif,
            args.flank,
            outdir / f"{args.output_prefix}.motif_centered.{args.flank}bp_flank.bed",
        )
        motif_center = motif_region.center
        region = read_first_bed_region(bed)
        LOGGER.info("Using motif-centered BED: %s", bed)
        LOGGER.info("Motif center: %s:%s", motif_region.chrom, motif_center)
        LOGGER.info("Plot region: %s:%s-%s", plot_region.chrom, plot_region.start, plot_region.end)
    else:
        bed = Path(args.bed).expanduser()
        if not bed.exists():
            raise FileNotFoundError(f"BED file not found: {bed}")
        region = read_first_bed_region(bed)

    center = args.cluster_center
    if center is None:
        center = motif_center if motif_center is not None else region.center

    samples = collect_samples(args)
    samples_by_condition = group_samples_by_condition(samples, args.condition)
    motif_markers = parse_marker_arguments(args.mark_motif)
    import_analysis_dependencies(require_scipy=args.cluster_waterfall == "hierarchical")

    if motif_markers:
        LOGGER.info("Motif markers:")
        for marker in motif_markers:
            LOGGER.info(
                "  %s: %s:%s-%s center=%s color=%s",
                marker.name,
                marker.chrom,
                marker.start,
                marker.end,
                marker.center,
                marker.color,
            )

    pooled_by_condition: "OrderedDict[str, pd.DataFrame]" = OrderedDict()

    for condition, condition_samples in samples_by_condition.items():
        LOGGER.info("Condition: %s", condition)
        matrices = []
        sample_names = []

        for sample in condition_samples:
            LOGGER.info("  Loading %s: %s", sample.sample, sample.bam)
            df = extract_matrix(
                bam=sample.bam,
                ref=ref,
                bed=bed,
                region_index=args.region_index,
                padding=args.padding,
                min_reads=args.min_reads,
            )
            matrices.append(df)
            sample_names.append(sample.sample)

        pooled = pd.concat(matrices, axis=0, keys=sample_names, names=["sample", "read"])
        pooled = pooled.reset_index(level=0, drop=True)
        pooled = sort_matrix_columns(pooled)
        pooled_by_condition[condition] = pooled

        matrix_out = outdir / f"{args.output_prefix}.{condition}.pooled_matrix.tsv"
        pooled.to_csv(matrix_out, sep="\t", index=True)
        LOGGER.info("  Pooled reads: %s", pooled.shape[0])
        LOGGER.info("  Matrix columns: %s", pooled.shape[1])
        LOGGER.info("  Matrix: %s", matrix_out)

    if args.cluster_waterfall != "none":
        LOGGER.info(
            "Waterfall clustering: %s, center=%s, window=+/- %s bp",
            args.cluster_waterfall,
            center,
            args.cluster_window,
        )

    waterfall_out = outdir / f"{args.output_prefix}.stacked_waterfalls_by_condition.{args.waterfall_format}"
    plot_stacked_waterfalls(
        condition_matrices=pooled_by_condition,
        out_path=waterfall_out,
        title=f"Pooled Fiber-seq waterfalls by condition\n{region.label}",
        center=center,
        cluster_window=args.cluster_window,
        cluster_method=args.cluster_waterfall,
        max_reads=args.max_waterfall_reads,
        motif_markers=motif_markers,
        random_seed=args.random_seed,
        dpi=args.dpi,
    )

    metaplot_out = outdir / f"{args.output_prefix}.metaplot_by_condition.{args.metaplot_format}"
    plot_condition_metaplot(
        condition_matrices=pooled_by_condition,
        out_path=metaplot_out,
        title=f"Pooled average profiles by condition\n{region.label}",
        smooth_window=args.smooth_window,
        motif_markers=motif_markers,
        ylabel=args.metaplot_ylabel,
        dpi=args.dpi,
    )

    summary_out = outdir / f"{args.output_prefix}.run_summary.txt"
    write_run_summary(summary_out, args, region, center, samples_by_condition)

    LOGGER.info("Done.")
    LOGGER.info("Output directory: %s", outdir)
    LOGGER.info("Stacked waterfall: %s", waterfall_out)
    LOGGER.info("Metaplot: %s", metaplot_out)
    LOGGER.info("Run summary: %s", summary_out)

    return {
        "outdir": outdir,
        "waterfall": waterfall_out,
        "metaplot": metaplot_out,
        "summary": summary_out,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose)

    try:
        run(args)
    except Exception as exc:
        LOGGER.error("%s", exc)
        if args.verbose:
            raise
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

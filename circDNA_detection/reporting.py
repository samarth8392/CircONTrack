#!/usr/bin/env python3
"""
Reporting and plotting helpers for CircONTrack candidate outputs.

The functions in this module accept structured candidates and simple coverage
window dictionaries. They do not rerun detection or mutate scientific scores.
"""

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from .confidence_scorer import ConfidenceScorer
from .utils import CircularCandidate


OKABE_ITO = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
]


def candidate_to_record(candidate: CircularCandidate, rank: int) -> Dict[str, Any]:
    """Convert a candidate object to a stable report/plot record."""
    return {
        "rank": rank,
        "chromosome": candidate.chromosome,
        "start": candidate.start,
        "end": candidate.end,
        "length": candidate.length,
        "confidence_score": candidate.confidence_score,
        "detection_method": candidate.detection_method,
        "mean_coverage": candidate.mean_coverage,
        "fold_enrichment": candidate.fold_enrichment,
        "coverage_uniformity": candidate.coverage_uniformity,
        "junction_support": candidate.junction_support,
        "split_support": candidate.split_support,
        "gc_content": candidate.gc_content,
    }


def candidates_to_records(candidates: Sequence[CircularCandidate]) -> List[Dict[str, Any]]:
    """Convert candidates to deterministic dictionaries sorted by rank."""
    return [candidate_to_record(candidate, i + 1) for i, candidate in enumerate(candidates)]


def confidence_breakdown(
    candidate: CircularCandidate,
    scorer: Optional[ConfidenceScorer] = None,
) -> Dict[str, float]:
    """Return score components using the same normalization as ConfidenceScorer."""
    scorer = scorer or ConfidenceScorer()
    components: Dict[str, float] = {}

    if candidate.fold_enrichment is not None:
        components["coverage_fold_enrichment"] = (
            min(candidate.fold_enrichment / scorer.thresholds["max_fold_enrichment"], 1.0)
            * scorer.weights["coverage"]["fold_enrichment"]
        )
    if candidate.coverage_uniformity is not None:
        components["coverage_uniformity"] = (
            max(0.0, candidate.coverage_uniformity)
            * scorer.weights["coverage"]["coverage_uniformity"]
        )
    if candidate.mean_coverage is not None:
        components["coverage_depth"] = (
            min(candidate.mean_coverage / scorer.thresholds["max_coverage"], 1.0)
            * scorer.weights["coverage"]["mean_coverage"]
        )
    if candidate.junction_support is not None:
        components["junction_support"] = (
            min(candidate.junction_support / scorer.thresholds["max_junction_support"], 1.0)
            * scorer.weights["junction"]["junction_support"]
        )
        components["junction_base"] = scorer.weights["junction"]["base_confidence"]
    if candidate.split_support is not None:
        components["split_support"] = (
            min(candidate.split_support / scorer.thresholds["max_split_support"], 1.0)
            * scorer.weights["split_read"]["split_support"]
        )
        components["split_base"] = scorer.weights["split_read"]["base_confidence"]

    components["final_confidence"] = candidate.confidence_score
    return components


def write_markdown_report(
    candidates: Sequence[CircularCandidate],
    output_file: str,
    run_metadata: Optional[Mapping[str, Any]] = None,
) -> Path:
    """Write a lightweight Markdown report for final CircONTrack candidates."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    records = candidates_to_records(candidates)
    metadata = dict(run_metadata or {})

    with output_path.open("w") as handle:
        handle.write("# CircONTrack Candidate Report\n\n")
        handle.write("## Overview\n\n")
        handle.write(
            "This report summarizes circular DNA candidates emitted by the implemented "
            "CircONTrack pipeline. It reports computed evidence fields and does not "
            "add benchmarked accuracy claims.\n\n"
        )

        if metadata:
            handle.write("## Run Configuration\n\n")
            for key in sorted(metadata):
                handle.write(f"- {key}: {metadata[key]}\n")
            handle.write("\n")

        handle.write("## Candidate Summary\n\n")
        handle.write(f"Total candidates after filtering: {len(records)}\n\n")

        handle.write("## Evidence Table\n\n")
        handle.write(
            "| Rank | Coordinate | Length (bp) | Method | Confidence | "
            "Fold enrichment | Junction reads | Split reads | GC content |\n"
        )
        handle.write("|---:|---|---:|---|---:|---:|---:|---:|---:|\n")
        if records:
            for record in records:
                coordinate = f"{record['chromosome']}:{record['start']}-{record['end']}"
                fold = _format_optional(record["fold_enrichment"], ".3f")
                junction = _format_optional(record["junction_support"], "d")
                split = _format_optional(record["split_support"], "d")
                gc = _format_optional(record["gc_content"], ".3f")
                handle.write(
                    f"| {record['rank']} | {coordinate} | {record['length']} | "
                    f"{record['detection_method']} | {record['confidence_score']:.3f} | "
                    f"{fold} | {junction} | {split} | {gc} |\n"
                )
        else:
            handle.write("| NA | No candidates passed filtering | NA | NA | NA | NA | NA | NA | NA |\n")
        handle.write("\n")

        handle.write("## Confidence Score Breakdown\n\n")
        if records:
            scorer = ConfidenceScorer()
            for record, candidate in zip(records, candidates):
                handle.write(
                    f"### Candidate {record['rank']}: "
                    f"{record['chromosome']}:{record['start']}-{record['end']}\n\n"
                )
                for key, value in confidence_breakdown(candidate, scorer).items():
                    handle.write(f"- {key}: {value:.3f}\n")
                handle.write("\n")
        else:
            handle.write("No confidence components were calculated because no candidates passed filtering.\n\n")

        handle.write("## Coordinate Convention\n\n")
        handle.write(
            "Candidate intervals and BED output use 0-based, half-open coordinates "
            "`[start, end)`, matching pysam fetch and BED conventions.\n\n"
        )

        handle.write("## Plot Interpretation\n\n")
        handle.write(
            "Summary plots visualize the final candidate table only. They are intended "
            "for inspection of evidence composition and ranking, not as independent "
            "statistical validation.\n"
        )

    return output_path


def plot_candidate_summary(
    candidates: Sequence[CircularCandidate],
    output_file: str,
    title: str = "CircONTrack Candidate Confidence Summary",
):
    """Plot final candidate confidence scores."""
    plt = _load_pyplot()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    _apply_plot_style(ax)

    records = candidates_to_records(candidates)
    if records:
        labels = [f"{r['chromosome']}:{r['start']}-{r['end']}" for r in records]
        scores = [r["confidence_score"] for r in records]
        ax.bar(range(len(records)), scores, color=OKABE_ITO[0])
        ax.set_xticks(range(len(records)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel("Confidence score (0-1)")
    else:
        _plot_empty_message(ax, "No candidates passed filtering")
        ax.set_ylabel("Confidence score (0-1)")

    ax.set_title(title)
    ax.set_xlabel("Candidate interval")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_evidence_support(
    candidates: Sequence[CircularCandidate],
    output_file: str,
    title: str = "CircONTrack Candidate Evidence Support",
):
    """Plot junction and split-read support for final candidates."""
    plt = _load_pyplot()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    _apply_plot_style(ax)

    records = candidates_to_records(candidates)
    if records:
        labels = [f"{r['chromosome']}:{r['start']}-{r['end']}" for r in records]
        x = list(range(len(records)))
        junction = [r["junction_support"] or 0 for r in records]
        split = [r["split_support"] or 0 for r in records]
        width = 0.38
        ax.bar([i - width / 2 for i in x], junction, width=width, label="Junction reads", color=OKABE_ITO[1])
        ax.bar([i + width / 2 for i in x], split, width=width, label="Split reads", color=OKABE_ITO[2])
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.legend(frameon=False)
    else:
        _plot_empty_message(ax, "No candidates passed filtering")

    ax.set_title(title)
    ax.set_xlabel("Candidate interval")
    ax.set_ylabel("Supporting read count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_confidence_breakdown(
    candidates: Sequence[CircularCandidate],
    output_file: str,
    title: str = "CircONTrack Confidence Score Components",
):
    """Plot score-component contributions for final candidates."""
    plt = _load_pyplot()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    _apply_plot_style(ax)

    records = candidates_to_records(candidates)
    component_names = [
        "coverage_fold_enrichment",
        "coverage_uniformity",
        "coverage_depth",
        "junction_support",
        "junction_base",
        "split_support",
        "split_base",
    ]

    if records:
        scorer = ConfidenceScorer()
        labels = [f"{r['chromosome']}:{r['start']}-{r['end']}" for r in records]
        bottoms = [0.0] * len(records)
        for i, component in enumerate(component_names):
            values = [confidence_breakdown(c, scorer).get(component, 0.0) for c in candidates]
            ax.bar(
                range(len(records)),
                values,
                bottom=bottoms,
                label=component.replace("_", " "),
                color=OKABE_ITO[i % len(OKABE_ITO)],
            )
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        ax.set_xticks(range(len(records)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, max(1.0, max(bottoms) * 1.1))
        ax.legend(frameon=False, fontsize=8)
    else:
        _plot_empty_message(ax, "No candidates passed filtering")

    ax.set_title(title)
    ax.set_xlabel("Candidate interval")
    ax.set_ylabel("Score contribution")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_coverage_profile(
    coverage_windows: Sequence[Mapping[str, Any]],
    candidates: Sequence[CircularCandidate],
    chromosome: str,
    output_file: str,
):
    """Plot a coverage profile from precomputed windows and candidate intervals."""
    plt = _load_pyplot()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    _apply_plot_style(ax)

    if coverage_windows:
        starts = [int(window["start"]) for window in coverage_windows]
        coverage = [float(window["coverage"]) for window in coverage_windows]
        ax.plot(starts, coverage, color=OKABE_ITO[0], linewidth=1.2, label="Window coverage")
        for candidate in candidates:
            if candidate.chromosome == chromosome:
                ax.axvspan(candidate.start, candidate.end, alpha=0.25, color=OKABE_ITO[1])
        ax.legend(frameon=False)
    else:
        _plot_empty_message(ax, "No coverage windows available")

    ax.set_title(f"Coverage Profile: {chromosome}")
    ax.set_xlabel(f"Position on {chromosome} (bp)")
    ax.set_ylabel("Window-normalized read count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def write_plot_bundle(
    candidates: Sequence[CircularCandidate],
    output_dir: str,
) -> List[Path]:
    """Write standard candidate summary plots and return created paths."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        out_dir / "candidate_confidence.png",
        out_dir / "candidate_evidence_support.png",
        out_dir / "confidence_breakdown.png",
    ]
    plot_candidate_summary(candidates, str(paths[0]))
    plot_evidence_support(candidates, str(paths[1]))
    plot_confidence_breakdown(candidates, str(paths[2]))
    return paths


def _format_optional(value: Any, fmt: str) -> str:
    if value is None:
        return "NA"
    if fmt == "d":
        return f"{int(value):d}"
    return format(float(value), fmt)


def _load_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "axes.prop_cycle": plt.cycler(color=OKABE_ITO),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })
    return plt


def _apply_plot_style(ax):
    ax.grid(True, axis="y", alpha=0.25)
    ax.set_axisbelow(True)


def _plot_empty_message(ax, message: str):
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])

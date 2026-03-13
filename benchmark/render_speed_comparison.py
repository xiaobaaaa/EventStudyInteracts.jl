from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

ENGINE_ORDER = [
    "Julia EventStudyInteracts.jl",
    "R fixest",
    "Stata eventstudyinteract",
]
COLORS = {
    "Julia EventStudyInteracts.jl": "#1b9e77",
    "R fixest": "#d95f02",
    "Stata eventstudyinteract": "#4c78a8",
}


def parse_cli(argv: list[str]) -> dict[str, str]:
    options: dict[str, str] = {}
    for arg in argv:
        if not arg.startswith("--") or "=" not in arg:
            continue
        key, value = arg[2:].split("=", 1)
        options[key] = value
    return options


def read_results(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def merge_results(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        rows.extend(read_results(path))
    return rows


def write_latest_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = ["engine", "spec", "seconds", "nobs"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


def escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def render_svg(path: Path, rows: list[dict[str, str]]) -> None:
    by_spec: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        row = dict(row)
        row["seconds"] = float(row["seconds"])
        by_spec[row["spec"]].append(row)

    specs = ["id + t", "id + id1 + id2 + t"]
    max_seconds = max(float(row["seconds"]) for row in rows)

    width = 1080
    height = 430
    margin_left = 240
    margin_right = 80
    margin_top = 88
    plot_width = width - margin_left - margin_right
    bar_height = 22
    bar_gap = 12
    group_gap = 52
    scale = plot_width / max_seconds if max_seconds else 1.0

    parts: list[str] = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append('<text x="60" y="42" font-size="26" font-family="Segoe UI, Arial, sans-serif" font-weight="700" fill="#1f2937">Event Study Speed Comparison</text>')
    parts.append('<text x="60" y="67" font-size="14" font-family="Segoe UI, Arial, sans-serif" fill="#4b5563">Synthetic Sun-Abraham benchmark on the same panel. Lower is better; Julia excludes first-run JIT compilation.</text>')

    legend_x = 60
    legend_y = 96
    for idx, engine in enumerate(ENGINE_ORDER):
        y = legend_y + idx * 24
        parts.append(f'<rect x="{legend_x}" y="{y - 11}" width="14" height="14" rx="3" fill="{COLORS[engine]}"/>')
        parts.append(f'<text x="{legend_x + 24}" y="{y}" font-size="13" font-family="Segoe UI, Arial, sans-serif" fill="#374151">{escape(engine)}</text>')

    baseline_y = height - 42
    parts.append(f'<line x1="{margin_left}" y1="{baseline_y}" x2="{width - margin_right}" y2="{baseline_y}" stroke="#d1d5db" stroke-width="1"/>')
    for tick in range(0, int(math.ceil(max_seconds / 10.0)) * 10 + 1, 10):
        x = margin_left + tick * scale
        parts.append(f'<line x1="{x:.1f}" y1="{margin_top - 12}" x2="{x:.1f}" y2="{baseline_y}" stroke="#f3f4f6" stroke-width="1"/>')
        parts.append(f'<text x="{x:.1f}" y="{baseline_y + 20}" text-anchor="middle" font-size="12" font-family="Segoe UI, Arial, sans-serif" fill="#6b7280">{tick}s</text>')

    current_y = 160
    for spec in specs:
        parts.append(f'<text x="60" y="{current_y + 16}" font-size="18" font-family="Segoe UI, Arial, sans-serif" font-weight="700" fill="#111827">{escape(spec)}</text>')
        rows_for_spec = {row["engine"]: row for row in by_spec[spec]}
        stata_seconds = rows_for_spec.get("Stata eventstudyinteract", {"seconds": 0.0})["seconds"]
        fastest = min(row["seconds"] for row in rows_for_spec.values())
        best_label = next(engine for engine, row in rows_for_spec.items() if row["seconds"] == fastest)
        parts.append(f'<text x="60" y="{current_y + 36}" font-size="12" font-family="Segoe UI, Arial, sans-serif" fill="#6b7280">Fastest: {escape(best_label)} | Julia speedup vs Stata: {stata_seconds / rows_for_spec["Julia EventStudyInteracts.jl"]["seconds"]:.2f}x</text>')

        bar_start_y = current_y + 52
        for idx, engine in enumerate(ENGINE_ORDER):
            row = rows_for_spec[engine]
            y = bar_start_y + idx * (bar_height + bar_gap)
            width_px = row["seconds"] * scale
            parts.append(f'<text x="{margin_left - 16}" y="{y + 16}" text-anchor="end" font-size="13" font-family="Segoe UI, Arial, sans-serif" fill="#374151">{escape(engine.replace(" EventStudyInteracts.jl", ""))}</text>')
            parts.append(f'<rect x="{margin_left}" y="{y}" width="{width_px:.1f}" height="{bar_height}" rx="6" fill="{COLORS[engine]}"/>')
            parts.append(f'<text x="{margin_left + width_px + 8:.1f}" y="{y + 16}" font-size="12" font-family="Segoe UI, Arial, sans-serif" fill="#111827">{row["seconds"]:.2f}s</text>')
        current_y += 3 * (bar_height + bar_gap) + group_gap

    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def render_summary(path: Path, rows: list[dict[str, str]]) -> None:
    by_spec: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        by_spec[row["spec"]][row["engine"]] = float(row["seconds"])

    lines = [
        "# Benchmark Summary",
        "",
        "Generated from the local Julia / R / Stata benchmark workflow.",
        "",
    ]
    for spec in ["id + t", "id + id1 + id2 + t"]:
        spec_rows = by_spec[spec]
        lines.append(f"## {spec}")
        lines.append("")
        for engine in ENGINE_ORDER:
            lines.append(f"- {engine}: {spec_rows[engine]:.3f} seconds")
        lines.append(f"- Julia vs Stata speedup: {spec_rows['Stata eventstudyinteract'] / spec_rows['Julia EventStudyInteracts.jl']:.2f}x")
        lines.append(f"- Julia vs fixest speedup: {spec_rows['R fixest'] / spec_rows['Julia EventStudyInteracts.jl']:.2f}x")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    import sys

    options = parse_cli(sys.argv[1:])
    input_paths = [Path(part) for part in options.get("inputs", "").split(",") if part]
    if not input_paths:
        raise SystemExit("Missing --inputs=...")

    svg_path = Path(options.get("svg", "benchmark/speed_comparison.svg"))
    csv_path = Path(options.get("csv", "benchmark/results_latest.csv"))
    summary_path = Path(options.get("summary", "benchmark/results_latest.md"))

    rows = merge_results(input_paths)
    write_latest_csv(csv_path, rows)
    render_svg(svg_path, rows)
    render_summary(summary_path, rows)
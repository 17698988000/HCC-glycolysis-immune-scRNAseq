from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.pdfgen import canvas


DARK = colors.HexColor("#243447")
GREY = colors.HexColor("#6B7280")
SAMPLE_COLORS = {
    "HCC1R": colors.HexColor("#4477AA"),
    "HCC2R": colors.HexColor("#EE6677"),
    "HCC3R": colors.HexColor("#228833"),
    "HCC4R": colors.HexColor("#CCBB44"),
}


def label(c, x, y, text, size=8, bold=False, color=DARK):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, text)


def draw_panel(c, data: pd.DataFrame, x: float, y: float, w: float, h: float, letter: str, title: str):
    label(c, x, y + h - 12, letter, 12, True)
    label(c, x + 15, y + h - 12, title, 8.5, True)
    outcomes = ["MIF", "SPP1", "MIF_SPP1_ligand_score", "four_gene_score"]
    outcome_labels = {
        "MIF": "MIF",
        "SPP1": "SPP1",
        "MIF_SPP1_ligand_score": "MIF/SPP1 score",
        "four_gene_score": "Four-gene score",
    }
    plot_x0, plot_x1 = x + 100, x + w - 25
    plot_y0, plot_y1 = y + 42, y + h - 40
    low, high = -0.2, 0.9

    def xpos(value):
        return plot_x0 + (value - low) / (high - low) * (plot_x1 - plot_x0)

    c.setStrokeColor(GREY)
    c.setDash(2, 2)
    c.line(xpos(0), plot_y0, xpos(0), plot_y1)
    c.setDash()
    samples = ["HCC1R", "HCC2R", "HCC3R", "HCC4R"]
    for i, outcome in enumerate(outcomes):
        yy = plot_y1 - 22 - i * 48
        label(c, x + 4, yy - 2, outcome_labels[outcome], 7)
        part = data[data["outcome"] == outcome].set_index("sample")
        for sample_index, sample in enumerate(samples):
            if sample not in part.index:
                continue
            value = float(part.loc[sample, "spearman_rho_block_means"])
            p = float(part.loc[sample, "empirical_p_positive"])
            c.setFillColor(SAMPLE_COLORS[sample])
            c.circle(xpos(value), yy + (sample_index - 1.5) * 5, 2.7, fill=1, stroke=0)
            if p < 0.05:
                label(c, xpos(value) + 4, yy + (sample_index - 1.5) * 5 - 2, "*", 6.5, True)
        positive = int((part["spearman_rho_block_means"] > 0).sum())
        significant = int((part["empirical_p_positive"] < 0.05).sum())
        label(c, plot_x1 - 70, yy + 12, f"{positive}/4 positive; {significant}/4 p<0.05", 6, color=GREY)
    c.setStrokeColor(GREY)
    c.line(plot_x0, plot_y0, plot_x1, plot_y0)
    for tick in [-0.2, 0, 0.4, 0.8]:
        xx = xpos(tick)
        c.line(xx, plot_y0 - 2, xx, plot_y0 + 2)
        label(c, xx - 8, plot_y0 - 13, f"{tick:.1f}", 6, color=GREY)
    label(c, plot_x0 + 70, plot_y0 - 25, "Spearman rho across spatial-block means", 6.5, color=GREY)
    for index, sample in enumerate(samples):
        xx = x + 20 + index * 65
        c.setFillColor(SAMPLE_COLORS[sample])
        c.circle(xx, y + 15, 2.5, fill=1, stroke=0)
        label(c, xx + 5, y + 12, sample, 6.5)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(args.results_csv)

    width, height = 900, 430
    c = canvas.Canvas(str(args.output), pagesize=(width, height))
    c.setFillColor(colors.white)
    c.rect(0, 0, width, height, fill=1, stroke=0)
    gap = 24
    panel_w = (width - 55 - gap) / 2
    panel_h = height - 65
    draw_panel(
        c,
        results[results["score"] == "locked_22_gene_mean_logexpr"],
        20,
        40,
        panel_w,
        panel_h,
        "A",
        "Locked 22-gene spatial score",
    )
    draw_panel(
        c,
        results[results["score"] == "leave_four_out_18_gene_mean_logexpr"],
        20 + panel_w + gap,
        40,
        panel_w,
        panel_h,
        "B",
        "Leave-four-out 18-gene spatial score",
    )
    label(
        c,
        20,
        12,
        "GSE238264 spatial-block empirical sensitivity: fixed 6x6 coordinate blocks, 10,000 within-section permutations; not RCTD-adjusted.",
        7,
        True,
    )
    c.showPage()
    c.save()
    print(args.output)


if __name__ == "__main__":
    main()

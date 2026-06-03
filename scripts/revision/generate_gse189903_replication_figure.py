from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.pdfgen import canvas


DARK = colors.HexColor("#243447")
GREY = colors.HexColor("#6B7280")
BLUE = colors.HexColor("#3B82B6")
RED = colors.HexColor("#C94C4C")
TEAL = colors.HexColor("#3C8D82")
PATIENT_COLORS = [
    colors.HexColor("#4477AA"),
    colors.HexColor("#EE6677"),
    colors.HexColor("#228833"),
    colors.HexColor("#CCBB44"),
]


def label(c, x, y, text, size=8, bold=False, color=DARK):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, text)


def draw_panel(c, data: pd.DataFrame, x: float, y: float, w: float, h: float, letter: str, title: str):
    label(c, x, y + h - 12, letter, 12, True)
    label(c, x + 15, y + h - 12, title, 9, True)
    features = ["MIF", "SPP1", "four_gene_readout"]
    feature_labels = {"MIF": "MIF", "SPP1": "SPP1", "four_gene_readout": "Four-gene readout"}
    plot_x0, plot_x1 = x + 110, x + w - 30
    plot_y0, plot_y1 = y + 38, y + h - 38
    all_values = data["effect_High_minus_Low"].astype(float)
    lim = max(abs(all_values.min()), abs(all_values.max()), 0.2) * 1.15

    def xpos(value):
        return plot_x0 + (value + lim) / (2 * lim) * (plot_x1 - plot_x0)

    c.setStrokeColor(GREY)
    c.setDash(2, 2)
    c.line(xpos(0), plot_y0, xpos(0), plot_y1)
    c.setDash()
    patients = sorted(data["patient"].unique())
    for i, feature in enumerate(features):
        yy = plot_y1 - 25 - i * 55
        label(c, x + 4, yy - 2, feature_labels[feature], 7.5)
        part = data[data["feature"] == feature].set_index("patient")
        for p_idx, patient in enumerate(patients):
            if patient not in part.index:
                continue
            value = float(part.loc[patient, "effect_High_minus_Low"])
            c.setFillColor(PATIENT_COLORS[p_idx % len(PATIENT_COLORS)])
            c.circle(xpos(value), yy + (p_idx - 1.5) * 4, 2.8, fill=1, stroke=0)
        positive = int((part["effect_High_minus_Low"] > 0).sum())
        label(c, plot_x1 - 28, yy + 10, f"{positive}/{len(part)} positive", 6.5, color=GREY)
    c.setStrokeColor(GREY)
    c.line(plot_x0, plot_y0, plot_x1, plot_y0)
    for tick in [-lim, 0, lim]:
        xx = xpos(tick)
        c.line(xx, plot_y0 - 2, xx, plot_y0 + 2)
        label(c, xx - 12, plot_y0 - 13, f"{tick:.2f}", 6, color=GREY)
    label(c, plot_x0 + 55, plot_y0 - 25, "Patient-level mean expression effect (GlycoHigh - GlycoLow)", 6.5, color=GREY)
    for p_idx, patient in enumerate(patients):
        xx = x + 15 + p_idx * 55
        c.setFillColor(PATIENT_COLORS[p_idx % len(PATIENT_COLORS)])
        c.circle(xx, y + 15, 2.5, fill=1, stroke=0)
        label(c, xx + 5, y + 12, patient, 6.5)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    effects = pd.read_csv(args.effects_csv)

    width, height = 900, 420
    c = canvas.Canvas(str(args.output), pagesize=(width, height))
    c.setFillColor(colors.white)
    c.rect(0, 0, width, height, fill=1, stroke=0)
    gap = 25
    panel_w = (width - 60 - gap) / 2
    panel_h = height - 65
    draw_panel(
        c,
        effects[effects["score"] == "locked_22_gene_rank_score"],
        25,
        40,
        panel_w,
        panel_h,
        "A",
        "Locked 22-gene rank score",
    )
    draw_panel(
        c,
        effects[effects["score"] == "leave_four_out_18_gene_rank_score"],
        25 + panel_w + gap,
        40,
        panel_w,
        panel_h,
        "B",
        "Leave-four-out 18-gene rank score",
    )
    label(
        c,
        25,
        12,
        "GSE189903 HCC-only independent single-cell replication: author-annotated malignant cells; no re-optimization.",
        7,
        True,
    )
    c.showPage()
    c.save()
    print(args.output)


if __name__ == "__main__":
    main()

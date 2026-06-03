from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.pdfgen import canvas


DARK = colors.HexColor("#243447")
GREY = colors.HexColor("#6B7280")
PATIENT_COLORS = [
    colors.HexColor("#4477AA"),
    colors.HexColor("#EE6677"),
    colors.HexColor("#228833"),
    colors.HexColor("#CCBB44"),
    colors.HexColor("#66CCEE"),
    colors.HexColor("#AA3377"),
    colors.HexColor("#BBBBBB"),
    colors.HexColor("#000000"),
]


def label(c, x, y, text, size=8, bold=False, color=DARK):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, text)


def draw_panel(
    c,
    data: pd.DataFrame,
    x: float,
    y: float,
    w: float,
    h: float,
    letter: str,
    title: str,
):
    label(c, x, y + h - 12, letter, 12, True)
    label(c, x + 15, y + h - 12, title, 8.2, True)
    features = ["MIF", "SPP1", "four_gene_readout"]
    feature_labels = {"MIF": "MIF", "SPP1": "SPP1", "four_gene_readout": "Four-gene readout"}
    plot_x0, plot_x1 = x + 92, x + w - 22
    plot_y0, plot_y1 = y + 34, y + h - 38
    all_values = data["effect_High_minus_Low"].astype(float)
    min_value = min(float(all_values.min()), 0.0)
    max_value = max(float(all_values.max()), 0.2)
    span = max_value - min_value
    low = min_value - span * 0.1
    high = max_value + span * 0.1

    def xpos(value):
        return plot_x0 + (value - low) / (high - low) * (plot_x1 - plot_x0)

    c.setStrokeColor(GREY)
    c.setDash(2, 2)
    c.line(xpos(0), plot_y0, xpos(0), plot_y1)
    c.setDash()
    patients = sorted(data["patient"].unique())
    for i, feature in enumerate(features):
        yy = plot_y1 - 24 - i * 55
        label(c, x + 4, yy - 2, feature_labels[feature], 7.2)
        part = data[data["feature"] == feature].set_index("patient")
        for p_idx, patient in enumerate(patients):
            if patient not in part.index:
                continue
            value = float(part.loc[patient, "effect_High_minus_Low"])
            c.setFillColor(PATIENT_COLORS[p_idx % len(PATIENT_COLORS)])
            c.circle(xpos(value), yy + (p_idx - 3.5) * 2.8, 2.2, fill=1, stroke=0)
        positive = int((part["effect_High_minus_Low"] > 0).sum())
        label(c, plot_x1 - 42, yy + 13, f"{positive}/{len(part)} positive", 6.2, color=GREY)
    c.setStrokeColor(GREY)
    c.line(plot_x0, plot_y0, plot_x1, plot_y0)
    for tick in [low, 0, high]:
        xx = xpos(tick)
        c.line(xx, plot_y0 - 2, xx, plot_y0 + 2)
        label(c, xx - 12, plot_y0 - 13, f"{tick:.2f}", 5.8, color=GREY)
    label(
        c,
        plot_x0 + 20,
        plot_y0 - 25,
        "Patient-level log-expression effect (GlycoHigh - GlycoLow)",
        6.2,
        color=GREY,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    effects = pd.read_csv(args.effects_csv)

    width, height = 1200, 420
    c = canvas.Canvas(str(args.output), pagesize=(width, height))
    c.setFillColor(colors.white)
    c.rect(0, 0, width, height, fill=1, stroke=0)
    gap = 18
    panel_w = (width - 50 - 2 * gap) / 3
    panel_h = height - 62
    panels = [
        ("locked_original_GlycoGroup", "A", "Original locked AUCell groups"),
        (
            "leave_four_out_global_rank_balanced",
            "B",
            "Leave-four-out global rank-balanced groups",
        ),
        (
            "leave_four_out_patient_internal_rank_balanced",
            "C",
            "Leave-four-out patient-internal groups",
        ),
    ]
    for index, (strategy, letter, title) in enumerate(panels):
        draw_panel(
            c,
            effects[effects["strategy"] == strategy],
            16 + index * (panel_w + gap),
            38,
            panel_w,
            panel_h,
            letter,
            title,
        )
    label(
        c,
        18,
        12,
        "GSE149614 discovery-cohort expression-level sensitivity: locked cells and genes; no coefficient or threshold optimization.",
        7,
        True,
    )
    c.showPage()
    c.save()
    print(args.output)


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape
from reportlab.pdfgen import canvas


BLUE = colors.HexColor("#3B82B6")
RED = colors.HexColor("#C94C4C")
TEAL = colors.HexColor("#3C8D82")
DARK = colors.HexColor("#243447")
GREY = colors.HexColor("#6B7280")
LIGHT = colors.HexColor("#E8EEF3")


def label(c: canvas.Canvas, x: float, y: float, text: str, size=8, bold=False, color=DARK):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, text)


def panel_title(c: canvas.Canvas, x: float, y: float, letter: str, title: str):
    label(c, x, y, letter, 12, True)
    label(c, x + 15, y, title, 9, True)


def coefficient_panel(c: canvas.Canvas, x: float, y: float, w: float, h: float):
    panel_title(c, x, y + h - 12, "A", "Four-gene tissue-level readout")
    genes = [("LDHA", 1.340), ("ENO1", 0.964), ("TPI1", 0.304), ("SLC2A1", 0.242)]
    max_val = 1.45
    base_y = y + h - 42
    for i, (gene, value) in enumerate(genes):
        yy = base_y - i * 24
        label(c, x + 4, yy + 3, gene, 8)
        bar_x = x + 55
        bar_w = (w - 90) * value / max_val
        c.setFillColor(TEAL)
        c.rect(bar_x, yy, bar_w, 10, fill=1, stroke=0)
        label(c, bar_x + bar_w + 4, yy + 2, f"{value:.3f}", 7)
    label(c, x + 4, y + 30, "Score = 0.304*TPI1 + 0.964*ENO1 + 1.340*LDHA + 0.242*SLC2A1", 6.5)
    label(c, x + 4, y + 17, "Discovery-stage tissue readout; not a clinical classifier.", 6.5, color=GREY)


def forest_panel(c: canvas.Canvas, x: float, y: float, w: float, h: float, letter: str, title: str, rows, note: str):
    panel_title(c, x, y + h - 12, letter, title)
    plot_x0 = x + 120
    plot_x1 = x + w - 76
    min_hr, max_hr = 0.5, 16.0

    def xpos(hr):
        import math

        return plot_x0 + (math.log(hr) - math.log(min_hr)) / (math.log(max_hr) - math.log(min_hr)) * (plot_x1 - plot_x0)

    axis_y = y + 35
    c.setStrokeColor(GREY)
    c.line(plot_x0, axis_y, plot_x1, axis_y)
    for tick in [0.5, 1, 2, 4, 8, 16]:
        xx = xpos(tick)
        c.line(xx, axis_y - 2, xx, axis_y + 2)
        label(c, xx - 4, axis_y - 12, str(tick), 6, color=GREY)
    c.setDash(2, 2)
    c.line(xpos(1), axis_y + 2, xpos(1), y + h - 30)
    c.setDash()
    label(c, plot_x0 + 18, axis_y - 23, "Hazard ratio (95% CI)", 6.5, color=GREY)

    base_y = y + h - 43
    for i, row in enumerate(rows):
        term, hr, lo, hi, p = row
        yy = base_y - i * 22
        label(c, x + 4, yy + 2, term, 7)
        c.setStrokeColor(DARK)
        c.setLineWidth(1)
        c.line(xpos(lo), yy + 4, xpos(hi), yy + 4)
        c.setFillColor(RED if term.startswith("Four-gene") else BLUE)
        c.circle(xpos(hr), yy + 4, 3, fill=1, stroke=0)
        label(c, plot_x1 + 7, yy + 1, f"{hr:.2f} ({lo:.2f}-{hi:.2f})", 6.5)
        label(c, plot_x1 + 61, yy + 1, p, 6.5)
    label(c, plot_x1 + 7, base_y + 18, "HR (95% CI)", 6.5, True)
    label(c, plot_x1 + 61, base_y + 18, "P", 6.5, True)
    label(c, x + 4, y + 10, note, 6.2, color=GREY)


def auc_panel(c: canvas.Canvas, x: float, y: float, w: float, h: float):
    panel_title(c, x, y + h - 12, "D", "TCGA-LIHC time-dependent AUC")
    values = [("1 year", 0.736), ("3 years", 0.701), ("5 years", 0.695)]
    chart_x = x + 35
    chart_y = y + 35
    chart_h = h - 70
    c.setStrokeColor(GREY)
    c.line(chart_x, chart_y, chart_x, chart_y + chart_h)
    c.line(chart_x, chart_y, x + w - 20, chart_y)
    for i, (term, value) in enumerate(values):
        bx = chart_x + 35 + i * 55
        bh = chart_h * value
        c.setFillColor(BLUE)
        c.rect(bx, chart_y, 28, bh, fill=1, stroke=0)
        label(c, bx + 1, chart_y + bh + 4, f"{value:.3f}", 7, True)
        label(c, bx - 2, chart_y - 13, term, 6.5)
    label(c, x + 4, y + 10, "Retrospective discrimination summary; not clinical validation.", 6.2, color=GREY)


def correlation_panel(c: canvas.Canvas, x: float, y: float, w: float, h: float):
    panel_title(c, x, y + h - 12, "E", "TCGA biological coverage")
    rows = [
        ("Core glycolysis", 0.76),
        ("SPP1/MIF axis", 0.49),
        ("Myeloid suppression", 0.46),
        ("T-cell dysfunction", 0.29),
        ("Hypoxia/angiogenesis", 0.22),
        ("Cytotoxic T-cell", 0.07),
    ]
    base_y = y + h - 38
    for i, (name, value) in enumerate(rows):
        yy = base_y - i * 17
        label(c, x + 4, yy + 2, name, 6.5)
        bar_x = x + 82
        bar_w = (w - 110) * value / 0.8
        c.setFillColor(TEAL)
        c.rect(bar_x, yy, bar_w, 8, fill=1, stroke=0)
        label(c, bar_x + bar_w + 3, yy + 1, f"{value:.2f}" + (" *" if value != 0.07 else " n.s."), 6.3)
    label(c, x + 4, y + 10, "* FDR < 0.05 across modules", 6.2, color=GREY)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    width, height = landscape((1104, 652))
    c = canvas.Canvas(str(args.output), pagesize=(width, height))
    c.setFillColor(colors.white)
    c.rect(0, 0, width, height, fill=1, stroke=0)

    margin = 25
    gap = 18
    top_h = 295
    bottom_h = 285
    col_w = (width - 2 * margin - 2 * gap) / 3

    coefficient_panel(c, margin, height - margin - top_h, col_w, top_h)
    forest_panel(
        c,
        margin + col_w + gap,
        height - margin - top_h,
        col_w,
        top_h,
        "B",
        "TCGA-LIHC multivariable Cox model",
        [
            ("Four-gene score, per SD", 1.59, 1.33, 1.91, "5.29e-7"),
            ("Age, per year", 1.01, 1.00, 1.03, "0.137"),
            ("Male vs female", 0.89, 0.61, 1.30, "0.537"),
            ("AJCC stage III-IV vs I-II", 1.81, 1.22, 2.68, "0.003"),
        ],
        "Clinical + score vs clinical only: LR p = 1.44e-06; C-index 0.704 vs 0.624.",
    )
    forest_panel(
        c,
        margin + 2 * (col_w + gap),
        height - margin - top_h,
        col_w,
        top_h,
        "C",
        "GSE14520 stable external evaluation",
        [
            ("Four-gene score, per SD", 1.39, 1.10, 1.74, "0.005"),
            ("AFP > 300 ng/mL", 1.54, 1.00, 2.38, "0.050"),
            ("Cirrhosis", 3.71, 0.90, 15.19, "0.069"),
            ("Multinodular disease", 1.29, 0.79, 2.11, "0.302"),
        ],
        "Stable binary AFP model: n = 218, events = 85; coefficients applied without re-optimization.",
    )
    auc_panel(c, margin + 80, margin, col_w + 40, bottom_h)
    correlation_panel(c, margin + col_w + gap + 110, margin, col_w + 90, bottom_h)
    label(
        c,
        margin,
        8,
        "Figure 6. State-derived four-gene tissue-level readout and retrospective cohort-level association.",
        7,
        True,
    )
    c.showPage()
    c.save()
    print(args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
benchmark_plots.py - benchmark figures for the PAPIO protocol, drawn from the tables
the README 9.4 benchmark already writes. No VCF parsing: every number in a figure is a
number in these tables, so text, tables and figures cannot disagree.

  somatic   somatic_metrics.tsv                             -> Figure3 (somatic accuracy, SEQC2)
  combined  germline_metrics.tsv + fusion_truth_result.tsv  -> Figure5 (germline HG002, fusions HCC1395)

usage:
  python3 benchmark_plots.py somatic  --somatic  $OUT_S/benchmark/somatic_metrics.tsv  [-o Figure3]
  python3 benchmark_plots.py combined --germline $OUT_G/benchmark/germline_metrics.tsv \
                                      --fusion   fusion_truth_result.tsv                [-o Figure5]

Metrics tables: columns callset type truth TP FP FN recall precision F1 [TP_call] (count_metrics output);
TP is the truth-side count; precision is computed from call-side counts (TP_call, FP).
Optional --final-dir <vcfeval_final_pass_main> adds panel F, SNV recall by truth VAF.
Output: <name>.pdf (vector), <name>.png (300 dpi), <name>.tiff (600 dpi). Needs pandas + matplotlib.
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

MM = 1 / 25.4
W = 174 * MM                                       # Cell Press full width
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "legend.frameon": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04})
# Okabe-Ito (colour-blind safe); every coloured mark is also labelled directly
BLUE, VERM, GREEN, PINK, ORANGE, SKY = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"
INK, INK2, GRID = "#222222", "#5A5A5A", "#E8E8E8"
CALLERS = [("mutect2_pass", "Mutect2", GREEN, "o"), ("strelka2_pass", "Strelka2", PINK, "s"),
           ("muse_pass", "MuSE", ORANGE, "^"), ("varscan_pass", "VarScan2", SKY, "v")]
STAGES = [("mutect2_raw", "Mutect2\nall calls"), ("mutect2_pass", "Filter-\nMutectCalls"),
          ("consensus_pass", "Consensus\n(>= 2 callers)"), ("final_pass", "Support\nfilter (final)")]


def die(msg):
    sys.exit(f"ERROR: {msg}")


def letter(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    for ext, dpi in (("pdf", 300), ("png", 300)):
        fig.savefig(f"{name}.{ext}", dpi=dpi)
    try:
        fig.savefig(f"{name}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception:
        pass
    plt.close(fig)
    print(f"wrote {os.path.abspath(name)}.pdf/.png/.tiff")


def read_metrics(path):
    if not os.path.isfile(path):
        die(f"file not found: {path}")
    t = pd.read_csv(path, sep=r"\s+")
    need = {"callset", "type", "truth", "TP", "FP", "FN", "recall", "precision", "F1"}
    if need - set(t.columns):
        die(f"{path} lacks columns {sorted(need - set(t.columns))}")
    # internal consistency: truth = TP + FN, recall = TP / truth (and precision = TP_call / (TP_call + FP))
    bad = t[(t.truth != t.TP + t.FN) | ((t.truth > 0) & ((t.recall - t.TP / t.truth.clip(lower=1)).abs() > 0.0006))]
    if "TP_call" in t.columns:
        n = (t.TP_call + t.FP).clip(lower=1)
        bad = pd.concat([bad, t[((t.TP_call + t.FP) > 0) & ((t.precision - t.TP_call / n).abs() > 0.0006)]])
    if len(bad):
        die(f"{path}: rows whose counts and rates disagree (re-run count_metrics):\n{bad.drop_duplicates().to_string(index=False)}")
    return t


def row(t, callset, vtype):
    r = t[(t.callset == callset) & (t.type == vtype)]
    if len(r) != 1:
        die(f"expected one row for {callset} / {vtype}, found {len(r)}")
    return r.iloc[0]


def vaf_values(path, key="TVAF"):
    """INFO/<key> of each SNV record (REF and ALT one base) in a vcfeval output VCF"""
    import gzip
    vals = []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            f = l.split("\t", 8)
            if len(f[3]) != 1 or any(len(x) != 1 for x in f[4].split(",")):
                continue
            m = re.search(rf"(?:^|;){key}=([^;]+)", f[7])
            if m:
                vals.append(float(m.group(1).split(",")[0]))
    return vals


# ======================================================================== somatic
def pr_panel(ax, t, vtype, title):
    """precision-recall of each caller's PASS calls and of the PAPIO steps (main scoring)"""
    pts = []
    for cs, name, col, mk in CALLERS:
        r = row(t, f"{cs}_main", vtype)
        if r.TP + r.FP == 0:                               # e.g. MuSE does not call indels
            continue
        pts.append((r.recall, r.precision, name, col, mk, 5.5, "caller"))
    for cs, name, mk, ms in (("consensus_pass_main", "Consensus", "D", 6), ("final_pass_main", "PAPIO final", "D", 7)):
        r = row(t, cs, vtype)
        pts.append((r.recall, r.precision, name, INK, mk, ms, "papio" if "final" in cs else "cons"))
    lo_x = max(0, np.floor((min(p[0] for p in pts) - 0.05) * 10) / 10)
    lo_y = max(0, np.floor((min(p[1] for p in pts) - 0.05) * 10) / 10)
    r_ = np.linspace(0.001, 1, 400)
    for f1 in (0.2, 0.4, 0.6, 0.8, 0.9):                   # iso-F1 curves
        p_ = f1 * r_ / (2 * r_ - f1)
        ok = (p_ > lo_y) & (p_ <= 1) & (r_ >= lo_x)
        if ok.sum() > 2 and r_[np.where(ok)[0][-1]] > 0.97:   # draw only curves that reach the right edge,
            ax.plot(r_[ok], p_[ok], color=GRID, lw=0.7, zorder=0)  # where they are labelled
            i = np.where(ok)[0][-1]
            if True:
                ax.text(r_[i], p_[i], f" F1 {f1:g}", fontsize=5.3, color="#A8A8A8", va="center", ha="left")
    for x, y, name, col, mk, ms, kind in pts:
        face = "white" if kind == "cons" else col
        ax.plot(x, y, mk, ms=ms, mfc=face, mec=col if kind != "cons" else INK, mew=1.0, zorder={"caller": 5, "cons": 4, "papio": 3}[kind],
                clip_on=False, label=f"{name}: {x + 1e-9:.3f} / {y + 1e-9:.3f}")
    ax.set_xlim(lo_x, 1.0); ax.set_ylim(lo_y, 1.0)
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_title(f"{title}  ({int(row(t, 'final_pass_main', vtype).truth):,} truth)", loc="left")
    order = {"PAPIO final": 0, "Consensus": 1, "Mutect2": 2, "Strelka2": 3, "MuSE": 4, "VarScan2": 5}
    h, l = ax.get_legend_handles_labels()
    hl = sorted(zip(h, l), key=lambda z: order[z[1].split(":")[0]])
    ax.legend([z[0] for z in hl], [z[1] for z in hl], loc="upper left", bbox_to_anchor=(-0.02, -0.24), ncol=2,
              title="Recall / precision", title_fontsize=6.5, alignment="left", columnspacing=1.0,
              handletextpad=0.3, fontsize=6.2)


def stage_panel(ax, t, vtype, title):
    """true and false calls left after each step (log scale)"""
    x = np.arange(len(STAGES))
    tp = [int(row(t, f"{cs}_main", vtype).TP) for cs, _ in STAGES]
    fp = [int(row(t, f"{cs}_main", vtype).FP) for cs, _ in STAGES]
    for y, col, mk, lab in ((tp, BLUE, "o", "Truth variants found"), (fp, VERM, "s", "False calls")):
        ax.plot(x, np.maximum(y, 0.8), "-", color=col, lw=1.6, zorder=2)
        ax.plot(x, np.maximum(y, 0.8), mk, ms=5, color=col, mec="white", mew=0.6, zorder=3, label=lab)
    for xi, a_, b_ in zip(x, tp, fp):
        for v, other in ((a_, b_), (b_, a_)):
            up = v >= other
            ax.text(xi, max(v, 0.8) * (1.6 if up else 0.6), f"{v:,}", ha="center",
                    va="bottom" if up else "top", fontsize=6, color=INK)
    ax.set_yscale("log"); ax.set_ylim(0.25, max(fp + tp) * 8)
    ax.set_xticks(x); ax.set_xticklabels([s for _, s in STAGES], fontsize=6)
    ax.set_xlim(-0.35, len(x) - 0.65)
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_ylabel("Count (log scale)")
    ax.set_title(title, loc="left")


def loss_panel(ax, t):
    """where the truth variants ended up: reported, or lost at which step"""
    cats = [("Reported (final PASS)", BLUE), ("Failed support filter", "#6E6E6E"),
            ("Failed consensus (< 2 callers)", "#9A9A9A"), ("Failed FilterMutectCalls", "#BDBDBD"),
            ("Not called by Mutect2", "#DEDEDE")]
    counts = []
    for i, (vt, lab) in enumerate((("snps", "SNVs"), ("indels", "Indels"))):
        raw, pas, con, fin = (int(row(t, f"{cs}_main", vt).TP) for cs, _ in STAGES)
        truth = int(row(t, "final_pass_main", vt).truth)
        vals = [fin, con - fin, pas - con, raw - pas, truth - raw]
        if sum(vals) != truth or min(vals) < 0:
            die(f"{vt}: stage counts are not nested (check the benchmark tables): {vals}, truth {truth}")
        counts.append(vals)
        left = 0
        y = 1 - i
        for v, (name, col) in zip(vals, cats):
            w = 100 * v / truth
            ax.barh(y, w, left=left, height=0.58, color=col, edgecolor="white", lw=1.0)
            if w >= 4:
                ax.text(left + w / 2, y, f"{v:,}", ha="center", va="center", fontsize=6,
                        color="white" if col in (BLUE, "#6E6E6E") else INK)
            left += w
    truths = [int(row(t, "final_pass_main", v).truth) for v in ("snps", "indels")]
    ax.set_yticks([1, 0]); ax.set_yticklabels([f"SNVs\n({truths[0]:,})", f"Indels\n({truths[1]:,})"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 100); ax.set_ylim(-0.5, 1.5); ax.set_xlabel("Truth variants (%)")
    ax.spines["left"].set_visible(False)
    ax.legend(handles=[Patch(color=c, label=f"{n}: {a_:,} / {b_:,}") for (n, c), a_, b_ in zip(cats, *counts)],
              ncol=2, loc="upper left", bbox_to_anchor=(-0.02, -0.32), handlelength=1.0, columnspacing=1.0, fontsize=6,
              title="SNVs / indels", title_fontsize=6.2, alignment="left")
    ax.set_title("Fate of the truth variants", loc="left")


def vaf_panel(ax, d):
    """SNV recall of the final call set by truth tumour VAF (INFO/TVAF of the SEQC2 truth records)"""
    tp, fn = (vaf_values(os.path.join(d, f)) for f in ("tp-baseline.vcf.gz", "fn.vcf.gz"))
    if not tp and not fn:
        die(f"no INFO/TVAF on SNV records in {d}/tp-baseline.vcf.gz or fn.vcf.gz")
    edges = [0, 0.05, 0.10, 0.20, 0.40, 1.0001]; labs = ["< 5", "5-10", "10-20", "20-40", ">= 40"]
    t, f = np.histogram(tp, bins=edges)[0], np.histogram(fn, bins=edges)[0]
    n = t + f
    rec = np.where(n > 0, t / np.maximum(n, 1), np.nan)
    x = np.arange(len(labs))
    ax.bar(x, 100 * rec, 0.62, color=BLUE, edgecolor="none")
    for xi, k, m, r in zip(x, t, n, rec):
        ax.text(xi, 100 * (r if r == r else 0) + 2, f"{100 * r:.0f}%\n{k}/{m}" if m else "0/0", ha="center",
                va="bottom", fontsize=5.6, color=INK2, linespacing=1.1)
    ax.set_xticks(x); ax.set_xticklabels(labs); ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 115); ax.set_yticks([0, 25, 50, 75, 100]); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel("Truth tumour VAF (%)"); ax.set_ylabel("Recall (%)")
    ax.set_title("Final SNV recall by allele fraction", loc="left")
    print("recall by VAF bin: " + ", ".join(f"{l}%: {k}/{m}" for l, k, m in zip(labs, t, n)))
    if sum(n) != int(row_cache["snv_truth"]):
        print(f"note: {sum(n)} SNVs with TVAF vs {row_cache['snv_truth']} truth SNVs in the table "
              "(records without TVAF are left out)")


row_cache = {}


def somatic(a):
    t = read_metrics(a.somatic)
    fig = plt.figure(figsize=(W, 175 * MM))
    # explicit layout (figure fractions) so the gaps between panel rows are even
    L, R, Wd = 0.08, 0.58, 0.39
    ax = fig.add_axes([L, 0.735, Wd, 0.235]); pr_panel(ax, t, "snps", "SNVs"); letter(ax, "A")
    ax = fig.add_axes([R, 0.735, Wd, 0.235]); pr_panel(ax, t, "indels", "Indels"); letter(ax, "B")
    ax = fig.add_axes([L, 0.385, Wd, 0.17]); stage_panel(ax, t, "snps", "SNVs after each step"); letter(ax, "C")
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.45))
    ax = fig.add_axes([R, 0.385, Wd, 0.17]); stage_panel(ax, t, "indels", "Indels after each step"); letter(ax, "D")
    if a.final_dir:
        ax = fig.add_axes([L, 0.115, Wd, 0.165]); loss_panel(ax, t); letter(ax, "E", x=-0.16, y=1.1)
        row_cache["snv_truth"] = row(t, "final_pass_main", "snps").truth
        ax = fig.add_axes([R, 0.115, Wd, 0.165]); vaf_panel(ax, a.final_dir); letter(ax, "F", x=-0.2, y=1.1)
    else:
        ax = fig.add_axes([L, 0.155, 0.89, 0.125]); loss_panel(ax, t); letter(ax, "E", x=-0.07, y=1.1)
    save(fig, a.o or "Figure3")
    f = row(t, "final_pass_main", "snps"); g = row(t, "final_pass_main", "indels")
    print(f"final SNV recall {f.recall:.3f} precision {f.precision:.3f} | indel recall {g.recall:.3f} "
          f"precision {g.precision:.3f} (strict {row(t, 'final_pass_strict', 'indels').precision:.3f})")


# ======================================================================== germline + fusion
CLONE = re.compile(r"(^|[-,(])(RP\d+-|RP\d+\.|AC\d{6}|AL\d{6}|AP\d{6}|CT[ABCD]-|LINC\d|LOC\d|KB-)", re.I)


def miss_reason(text):
    first = re.split(r",\s+", str(text).strip())[0]
    if re.fullmatch(r"([A-Za-z0-9.]+)-\1", first):
        return "Within one gene"
    if "(" in first:
        return "Intergenic breakpoint"
    if CLONE.search(first):
        return "Non-coding partner"
    return "Two coding genes"


def germline_panel(ax, t):
    rows = []
    for vt, lab in (("snps", "SNV"), ("indels", "Indel")):
        b, f = row(t, "germline_before_filter", vt), row(t, "germline_after_filter", vt)
        rows += [(f"{lab} recall", b.recall, f.recall), (f"{lab} precision", b.precision, f.precision)]
    y = np.arange(len(rows))[::-1]
    for yy, (lab, vb, va) in zip(y, rows):
        ax.annotate("", xy=(va, yy), xytext=(vb, yy),
                    arrowprops=dict(arrowstyle="-|>", color="#9A9A9A", lw=1.0, mutation_scale=7,
                                    shrinkA=3.5, shrinkB=3.5))
        ax.plot(vb, yy, "o", ms=5.5, mfc="white", mec=INK, mew=1.0, zorder=3)
        ax.plot(va, yy, "o", ms=5.5, color=INK, zorder=4)
    lo = np.floor((min(min(r[1], r[2]) for r in rows) - 0.03) * 20) / 20
    ax.set_xlim(max(0, lo), 1.0); ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(y); ax.set_yticklabels([f"{r[0]}\n{r[1] + 1e-9:.3f} ->{r[2] + 1e-9:.3f}" for r in rows], fontsize=6.3)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel("Recall or precision")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="Before filter"),
                       Line2D([], [], marker="o", ls="", color=INK, label="After filter (PASS)")],
              loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, columnspacing=1.0, handletextpad=0.3)
    truth = int(row(t, "germline_after_filter", "snps").truth), int(row(t, "germline_after_filter", "indels").truth)
    ax.set_title(f"Germline, HG002\n({truth[0]:,} SNVs, {truth[1]:,} indels)", loc="left")


def fusion_events(path):
    if not os.path.isfile(path):
        die(f"file not found: {path}")
    r = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    need = {"SV_type", "truth_fusion", "status", "matched_call", "TIER"}
    if need - set(r.columns):
        die(f"{path} lacks columns {sorted(need - set(r.columns))}")
    found = r[r.status == "FOUND"].drop_duplicates("matched_call")    # two reference rows can share one call
    missed = r[r.status != "FOUND"].copy()
    missed["reason"] = missed.truth_fusion.map(miss_reason)
    return found, missed, len(r)


def fusion_type_panel(ax, found, missed):
    order = ["TRA", "INV", "DUP", "INS", "DEL"]
    types = [s for s in order if s in set(found.SV_type) | set(missed.SV_type)]
    types += sorted((set(found.SV_type) | set(missed.SV_type)) - set(types))
    names = {"TRA": "Translocation", "INV": "Inversion", "DUP": "Duplication", "INS": "Insertion", "DEL": "Deletion"}
    y = np.arange(len(types))[::-1]
    for yy, s in zip(y, types):
        f, m = int((found.SV_type == s).sum()), int((missed.SV_type == s).sum())
        ax.barh(yy, f, height=0.6, color=BLUE, edgecolor="white", lw=1.0)
        ax.barh(yy, m, left=f, height=0.6, color="#D0D0D0", edgecolor="white", lw=1.0)
        ax.text(f + m + 0.25, yy, f"{f}/{f + m}", va="center", fontsize=6.2, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([names.get(s, s) for s in types]); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, max(((found.SV_type == s).sum() + (missed.SV_type == s).sum()) for s in types) + 2)
    ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel("Reference events")
    ax.legend(handles=[Patch(color=BLUE, label="Recovered"), Patch(color="#D0D0D0", label="Missed")],
              loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, handlelength=1.0)
    n = len(found) + len(missed)
    ax.set_title(f"Fusions, HCC1395\n{len(found)} of {n} recovered", loc="left")


def fusion_detail_panel(ax, found, missed):
    # FUSION_CONSENSUS verdicts: PASS_CONSENSUS (>= fusion_min_callers callers), PASS_KNOWN (1 caller,
    # known oncogenic pair), REVIEW (1 caller); anything else is shown under its own name
    names = {"PASS_CONSENSUS": ">= 2 callers", "PASS_KNOWN": "1 caller, known fusion", "REVIEW": "1 caller"}
    if "Verdict" not in found.columns:
        die("fusion result has no Verdict column")
    counts = found.Verdict.value_counts()
    order = [v for v in names if v in counts] + sorted(set(counts.index) - set(names))
    seg = [(names.get(v, v), int(counts[v])) for v in order]
    groups = []                                            # (row label, [(segment label, n)], colours)
    groups.append(("Recovered:\ncaller support", seg, [BLUE, SKY, "#A9CBE3", "#D6E6F2"][:len(seg)]))
    reasons = ["Two coding genes", "Non-coding partner", "Intergenic breakpoint", "Within one gene"]
    groups.append(("Missed:\nreference event", [(r, int((missed.reason == r).sum())) for r in reasons
                                                if (missed.reason == r).any()],
                   ["#5E5E5E", "#878787", "#ADADAD", "#D0D0D0"]))
    for yy, (lab, segs, cols) in zip((1, 0), groups):
        left = 0
        for (name, n), col in zip(segs, cols):
            ax.barh(yy, n, left=left, height=0.6, color=col, edgecolor="white", lw=1.0)
            ax.text(left + n / 2, yy, str(n), ha="center", va="center", fontsize=6.2,
                    color="white" if col in (BLUE, "#5E5E5E", "#878787") else INK)
            left += n
    ax.set_yticks([1, 0]); ax.set_yticklabels([g[0] for g in groups]); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, max(len(found), len(missed)) + 0.5); ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Reference events")
    hand = [Patch(color=c, label=n) for g in groups for (n, _), c in zip(g[1], g[2])]
    ax.legend(handles=hand, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, handlelength=1.0,
              columnspacing=0.8, fontsize=6)
    ax.set_title("Caller support and\nreasons for misses", loc="left")


def combined(a):
    t = read_metrics(a.germline)
    found, missed, nrows = fusion_events(a.fusion)
    fig = plt.figure(figsize=(W, 70 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1, 1], wspace=0.95)
    ax = fig.add_subplot(gs[0, 0]); germline_panel(ax, t); letter(ax, "A", x=-0.6)
    ax = fig.add_subplot(gs[0, 1]); fusion_type_panel(ax, found, missed); letter(ax, "B", x=-0.5)
    ax = fig.add_subplot(gs[0, 2]); fusion_detail_panel(ax, found, missed); letter(ax, "C", x=-0.5)
    save(fig, a.o or "Figure5")
    print(f"fusions: {len(found)} of {len(found) + len(missed)} distinct reference events recovered "
          f"({nrows} reference rows; {nrows - len(found) - len(missed)} share a call)")
    print("recovered by verdict:", found.Verdict.value_counts().to_dict())
    print("recovered by tier:", found.TIER.value_counts().to_dict())
    print("missed by reason:", missed.reason.value_counts().to_dict())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("somatic"); p.add_argument("--somatic", required=True); p.add_argument("-o")
    p.add_argument("--final-dir", help="vcfeval folder of the final call set (adds panel F, recall by truth VAF)")
    p = sp.add_parser("combined"); p.add_argument("--germline", required=True)
    p.add_argument("--fusion", required=True); p.add_argument("-o")
    a = ap.parse_args()
    {"somatic": somatic, "combined": combined}[a.cmd](a)


if __name__ == "__main__":
    main()
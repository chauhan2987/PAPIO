#!/usr/bin/env python3
"""
make_figures.py - PAPIO protocol figures drawn from pipeline output files.

  qc    Figure 2  sequencing and sample QC (tumor vs normal)
  cnv   Figure 4  copy number, arm-level gains/losses/LOH, purity and ploidy (Sequenza)

The benchmark figures (somatic accuracy, germline and fusions) are made by benchmark_plots.py,
and the resource table by trace_table.py.

You give the run folder and sample IDs; the script finds each file by its pipeline name in any
subfolder, and stops with a clear message if a file is missing or present more than once.

usage:
  python3 make_figures.py qc  --results $OUT_S --tumor valid1_T --normal valid1_N [-o Figure2]
  python3 make_figures.py cnv --results $OUT_S --tumor valid1_T                    [-o Figure4]

Files used
  qc : <id>.hs_metrics (Picard CollectHsMetrics), <id>.thresholds.bed.gz and <id>.mosdepth.summary.txt
       (mosdepth), <id>_dup_metrics.* (samtools markdup or Picard), somalier.pairs.tsv (optional),
       <tumor>.contamination.table (GATK CalculateContamination, optional)
  cnv: <tumor>_segments.txt and <tumor>_alternative_solutions.txt (Sequenza)

Output: <name>.pdf (vector), <name>.png (300 dpi), <name>.tiff (600 dpi). Needs pandas, numpy, matplotlib.
"""
import argparse, fnmatch, gzip, json, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MM = 1 / 25.4
W = 174 * MM                                       # Cell Press full width
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "legend.frameon": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04})
BLUE, VERM, GREY = "#0072B2", "#D55E00", "#8C8C8C"           # Okabe-Ito
SAMPLE = {"Tumor": VERM, "Normal": BLUE}
INK, INK2, GRID = "#222222", "#5A5A5A", "#E8E8E8"

# GRCh38: chromosome length and approximate centromere midpoint (inside the UCSC cytoBand acen band),
# used for genome coordinates and to split arms; a ~1 Mb shift does not change arm-level fractions materially
GRCH38 = {  # chrom: (length, centromere)
    "1": (248956422, 123400000), "2": (242193529, 93900000), "3": (198295559, 90900000),
    "4": (190214555, 50000000), "5": (181538259, 48800000), "6": (170805979, 59800000),
    "7": (159345973, 60100000), "8": (145138636, 45200000), "9": (138394717, 43000000),
    "10": (133797422, 39800000), "11": (135086622, 53400000), "12": (133275309, 35500000),
    "13": (114364328, 17700000), "14": (107043718, 17200000), "15": (101991189, 19000000),
    "16": (90338345, 36800000), "17": (83257441, 25100000), "18": (80373285, 18500000),
    "19": (58617616, 26200000), "20": (64444167, 28100000), "21": (46709983, 12000000),
    "22": (50818468, 15000000), "X": (156040895, 60600000), "Y": (57227415, 10400000)}
ACROCENTRIC = {"13", "14", "15", "21", "22"}       # p arms carry no unique sequence: not reported


def die(msg):
    sys.exit(f"\nERROR: {msg}\n")


def letter(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    for ext, dpi in (("pdf", 300), ("png", 300)):
        try:
            fig.savefig(f"{name}.{ext}", dpi=dpi)
        except OSError as e:
            die(f"cannot write {name}.{ext} in {os.getcwd()} ({e}); run from a folder you can write to")
    try:
        fig.savefig(f"{name}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception:
        pass
    plt.close(fig)
    print(f"wrote {os.path.abspath(name)}.pdf/.png/.tiff")


# ------------------------------------------------------------------ finding files
_files = {}


def find(root, name, required=True):
    """the one file called `name` (wildcards allowed) anywhere under root"""
    if root not in _files:
        if not os.path.isdir(root):
            die(f"folder not found: {root}\nIf it exists, the container cannot see it: bind it with -B {root}")
        _files[root] = [os.path.join(d, f) for d, _, fs in os.walk(root) for f in fs]
    hits = sorted(p for p in _files[root] if fnmatch.fnmatch(os.path.basename(p), name))
    if len(hits) == 1:
        return hits[0]
    if not hits:
        if required:
            die(f"no file named '{name}' under {root}")
        return None
    die(f"{len(hits)} files named '{name}' under {root}:\n  " + "\n  ".join(hits) +
        "\nGive --results the run folder that holds only one of them.")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_tsv(path):
    t = pd.read_csv(path, sep="\t", dtype=str)
    t.columns = [c.strip().strip('"').lstrip("#") for c in t.columns]
    return t.apply(lambda s: s.str.strip('"'))


# ================================================================== Figure 2: QC
def picard_rows(path):
    lines = [l.rstrip("\n") for l in opener(path)]
    for i, l in enumerate(lines):
        if l.startswith("## METRICS CLASS"):
            hdr, rows = lines[i + 1].split("\t"), []
            for x in lines[i + 2:]:
                if not x.strip() or x.startswith("#"):
                    break
                f = x.split("\t")
                rows.append(f + [""] * (len(hdr) - len(f)))
            return pd.DataFrame(rows, columns=hdr)
    die(f"{path} is not a Picard metrics file")


def hs_values(path):
    r = picard_rows(path).iloc[0]
    g = lambda k: float(r[k]) if k in r and r[k] not in ("", "?") else np.nan
    return {"mean": g("MEAN_TARGET_COVERAGE"), "median": g("MEDIAN_TARGET_COVERAGE"),
            "pct": {t: 100 * g(f"PCT_TARGET_BASES_{t}X") for t in (1, 10, 20, 30, 40, 50, 100)
                    if f"PCT_TARGET_BASES_{t}X" in r},
            "fold80": g("FOLD_80_BASE_PENALTY"), "on_target": 100 * g("PCT_SELECTED_BASES")}


def dup_rate(path):
    """fraction of primary reads marked duplicate: samtools markdup (text or --json) or Picard"""
    txt = opener(path).read()
    if "## METRICS CLASS" in txt:
        m = picard_rows(path)
        n = lambda c: pd.to_numeric(m[c], errors="coerce").sum()
        return (n("UNPAIRED_READ_DUPLICATES") + 2 * n("READ_PAIR_DUPLICATES")) / \
               (n("UNPAIRED_READS_EXAMINED") + 2 * n("READ_PAIRS_EXAMINED"))
    try:
        d = json.loads(txt)
    except ValueError:
        d = dict(l.split(":", 1) for l in txt.splitlines() if ":" in l)
    d = {k.strip(): float(str(v).strip()) for k, v in d.items() if str(v).strip().replace(".", "", 1).isdigit()}
    if d.get("EXAMINED", 0) > 0:
        return d.get("DUPLICATE PRIMARY TOTAL", d.get("DUPLICATE PAIR", 0) + d.get("DUPLICATE SINGLE", 0)) / d["EXAMINED"]
    die(f"{path}: not a samtools markdup or Picard duplicate report")


def target_breadth(path, depth=20):
    """% of each capture target covered at >= depth (mosdepth --thresholds)"""
    t = pd.read_csv(path, sep="\t")
    t.columns = [c.lstrip("#") for c in t.columns]
    if f"{depth}X" not in t.columns:
        die(f"{path} has no {depth}X column")
    return (100 * t[f"{depth}X"] / (t["end"] - t["start"]).clip(lower=1)).clip(0, 100).values


def sex_ratios(path):
    """chrX and chrY mean depth relative to autosomes, over the targets (_region rows)"""
    s = pd.read_csv(path, sep="\t")
    s["chrom"] = s["chrom"].astype(str)
    reg = s[s.chrom.str.endswith("_region")]
    if len(reg):
        s = reg.assign(chrom=reg.chrom.str.replace("_region$", "", regex=True))
    s["c"] = s.chrom.str.replace("^chr", "", regex=True)
    auto = s[s.c.isin([str(i) for i in range(1, 23)])]
    a = (auto["mean"] * auto["length"]).sum() / auto["length"].sum()
    get = lambda c: float(s.loc[s.c == c, "mean"].iloc[0]) / a if (s.c == c).any() else np.nan
    return get("X"), get("Y")


def qc(a):
    R, labs, ids = a.results, ["Tumor", "Normal"], [a.tumor, a.normal]
    hsf = [find(R, f"{i}.hs_metrics") for i in ids]
    thf = [find(R, f"{i}.thresholds.bed.gz") for i in ids]
    smf = [find(R, f"{i}.mosdepth.summary.txt") for i in ids]
    duf = [find(R, f"{i}_dup_metrics.*") for i in ids]
    pairs = find(R, "somalier.pairs.tsv", required=False)
    contam = find(R, f"{a.tumor}.contamination.table", required=False)
    hs, dups = [hs_values(p) for p in hsf], [dup_rate(p) for p in duf]
    br, sx = [target_breadth(p) for p in thf], [sex_ratios(p) for p in smf]

    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(W, 66 * MM), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    fig.subplots_adjust(wspace=0.85)

    # A - breadth: % of target bases covered at >= each depth (HsMetrics PCT_TARGET_BASES_nX), one bar per depth
    th = sorted(set(hs[0]["pct"]) & set(hs[1]["pct"]))
    xa, wa = np.arange(len(th)), 0.4
    for i, (lab, h, du) in enumerate(zip(labs, hs, dups)):
        a0.bar(xa + (i - 0.5) * wa, [h["pct"][t] for t in th], wa * 0.9, color=SAMPLE[lab], edgecolor="none",
               label=f"{lab}: mean {h['mean']:.1f}X, median {h['median']:.0f}X\n{100 * du:.2f}% duplicates")
    a0.set_xticks(xa); a0.set_xticklabels([f"{t}" for t in th]); a0.tick_params(axis="x", length=0)
    a0.set_ylim(0, 100); a0.yaxis.grid(True, color=GRID, lw=0.5); a0.set_axisbelow(True)
    a0.set_xlabel("Minimum depth (X)"); a0.set_ylabel("Target bases at >= depth (%)")
    a0.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.24), labelspacing=0.6, handlelength=1.0)
    a0.set_title("Breadth of target coverage", loc="left"); letter(a0, "A")


    # B - uniformity: share of capture targets by % of their bases at >= 20x (mosdepth)
    # y(x) = % of targets with at least x% of their bases at >= 20x  (exact, evaluated on a 0.1% grid)
    grid = np.round(np.arange(0, 100.05, 0.1), 1)
    for lab, b, h in zip(labs, br, hs):
        a1.plot(grid, [100 * np.mean(b >= g) for g in grid], color=SAMPLE[lab], lw=1.4,
                label=f"{lab}: {100 * np.mean(b >= 90):.2f}%\nfold-80 {h['fold80']:.2f}")
    a1.axvline(90, ls=(0, (3, 2)), lw=0.6, color=GREY, zorder=0)
    a1.set_xlim(0, 100); a1.set_ylim(0, 102); a1.yaxis.grid(True, color=GRID, lw=0.5); a1.set_axisbelow(True)
    a1.set_xlabel("Bases of the target at >= 20X (%)"); a1.set_ylabel("Capture targets (%)")
    a1.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.24), title="Targets >= 90% covered at 20X",
              title_fontsize=6.5, alignment="left", labelspacing=0.6, ncol=2, columnspacing=1.0, handlelength=1.2)
    a1.set_title(f"Uniformity across {len(br[0]):,} targets", loc="left"); letter(a1, "B")

    # C - sex check (chrX/chrY depth vs autosomes); identity and contamination against thresholds
    x, w = np.arange(2), 0.42
    for i, lab in enumerate(labs):
        a2.bar(x + (i - 0.5) * w, sx[i], w * 0.85, color=SAMPLE[lab], label=lab, edgecolor="none")
        for xx, v in zip(x + (i - 0.5) * w, sx[i]):
            if v == v:
                a2.text(xx, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=5.5, color=INK2)
    for yv, txt in ((1.0, "2 copies"), (0.5, "1 copy")):
        a2.axhline(yv, ls=(0, (3, 2)), lw=0.6, color=GREY, zorder=0)
        a2.text(1.62, yv, txt, fontsize=5.8, color=INK2, va="center", ha="left")
    a2.set_xticks(x); a2.set_xticklabels(["chrX", "chrY"]); a2.set_xlim(-0.55, 1.6)
    a2.set_ylabel("Depth relative to autosomes")
    notes, rel, cv = [], None, None
    if pairs:
        p = read_tsv(pairs)
        sel = p[((p.sample_a == a.tumor) & (p.sample_b == a.normal)) | ((p.sample_a == a.normal) & (p.sample_b == a.tumor))]
        if len(sel):

            # rel = float(sel.relatedness.iloc[0])
            # ok = "pass" if rel >= a.min_relatedness else "FAIL"
            # same rule as CHECK_SAMPLE_SWAPS: tumour LOH lowers relatedness (het -> hom) but never
            # creates opposite homozygotes, so IBS0/n ~ 0 means the same person despite low relatedness
            r0 = sel.iloc[0]
            rel = float(r0.relatedness)
            num = lambda k: pd.to_numeric(r0.get(k), errors="coerce")   # missing column / empty cell -> NaN
            f0 = num("ibs0") / num("n") if num("n") > 0 else np.nan
            if rel >= a.min_relatedness:
                ok = "pass"
            elif f0 == f0 and f0 <= a.max_ibs0:
                ok = f"IBS0 {f0:.4f}: same person, tumour LOH"
            else:
                ok = "FAIL"
            notes.append(f"Tumor-normal relatedness {rel:.2f} (>= {a.min_relatedness:g}: {ok})")
        else:
            print(f"note: pair {a.tumor}/{a.normal} not in {pairs}; relatedness left out")
    else:
        print("note: somalier.pairs.tsv not found; relatedness left out")
    if contam:
        cv = float(pd.read_csv(contam, sep="\t", comment="#")["contamination"].iloc[0])
        ok = "pass" if cv <= a.max_contamination else "FAIL"
        notes.append(f"Tumor contamination {100 * cv:.2f}% (<= {100 * a.max_contamination:g}%: {ok})")
    else:
        print(f"note: {a.tumor}.contamination.table not found; contamination left out")
    top = np.nanmax(np.array(sx, dtype=float))
    a2.set_ylim(0, max(1.55, top + 0.4))
    if notes:
        a2.text(0.04, 0.99, "\n".join(notes), transform=a2.transAxes, fontsize=5.8, va="top", ha="left",
                color=INK2, linespacing=1.5)
    a2.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.24), ncol=2, handlelength=1.0)
    a2.set_title("Sex, identity and contamination", loc="left"); letter(a2, "C")
    save(fig, a.o or "Figure2")

    print("\nValues in the figure (quote these in the text):")
    for lab, i, h, d, s in zip(labs, ids, hs, dups, sx):
        print(f"  {lab} ({i}): mean {h['mean']:.1f}x, median {h['median']:.1f}x, {h['pct'].get(30, np.nan):.1f}% "
              f"of target bases >= 30x, " + (f"on-target {h['on_target']:.1f}%, " if h['on_target'] == h['on_target'] else "") +
              f"duplicates {100 * d:.1f}%, "
              f"fold-80 {h['fold80']:.2f}, chrX/auto {s[0]:.2f}, chrY/auto {s[1]:.2f}")
    if rel is not None:
        print(f"  relatedness {rel:.3f}")
    if cv is not None:
        print(f"  contamination {100 * cv:.2f}%")
    print("Files used:\n  " + "\n  ".join(hsf + thf + smf + duf + [p for p in (pairs, contam) if p]))


# ================================================================== Figure 4: copy number
def arm_table(seg, base):
    """per autosomal arm: % of the segmented length that is gained, lost (vs base copy number) or LOH (B = 0)"""
    rows = []
    for c, (length, cen) in GRCH38.items():
        if c in ("X", "Y"):                        # sex chromosomes: expected copy number depends on sex
            continue
        s = seg[seg.chr == c]
        for arm, lo, hi in (("p", 0, cen), ("q", cen, length)):
            if arm == "p" and c in ACROCENTRIC:
                continue
            w = (s["end.pos"].clip(lo, hi) - s["start.pos"].clip(lo, hi)).clip(lower=0)
            if w.sum() == 0:
                continue
            rows.append(dict(arm=f"{c}{arm}", gain=100 * w[s.CNt > base].sum() / w.sum(),
                             loss=100 * w[s.CNt < base].sum() / w.sum(),
                             loh=100 * w[(s.B == 0) & (s.CNt > 0)].sum() / w.sum()))
    return pd.DataFrame(rows)


def cnv(a):
    segf = find(a.results, f"{a.tumor}_segments.txt")
    solf = find(a.results, f"{a.tumor}_alternative_solutions.txt")
    seg, sol = read_tsv(segf), read_tsv(solf)
    for c in ("chromosome", "start.pos", "end.pos", "CNt", "B"):
        if c not in seg.columns:
            die(f"{segf} has no '{c}' column")
    for c in ("cellularity", "ploidy", "SLPP"):
        if c not in sol.columns:
            die(f"{solf} has no '{c}' column")
    for c in ("start.pos", "end.pos", "CNt", "B"):
        seg[c] = pd.to_numeric(seg[c], errors="coerce")
    seg["chr"] = seg.chromosome.str.replace("^chr", "", regex=True)
    seg = seg[seg.chr.isin(GRCH38) & seg.CNt.notna()].copy()
    sol = sol[["cellularity", "ploidy", "SLPP"]].apply(pd.to_numeric, errors="coerce")
    sol = sol.sort_values("SLPP", ascending=False).reset_index(drop=True)
    cell, ploidy = float(sol.cellularity[0]), float(sol.ploidy[0])          # Sequenza's selected solution
    base = int(np.floor(ploidy + 0.5))                                    # round half up (2.5 -> 3)

    chroms = [c for c in GRCH38 if c in set(seg.chr)]
    lens = pd.Series({c: GRCH38[c][0] for c in chroms}, dtype=float)
    offs = lens.cumsum().shift(fill_value=0)
    seg["g0"] = seg["start.pos"] + seg.chr.map(offs)
    seg["g1"] = seg["end.pos"] + seg.chr.map(offs)

    fig = plt.figure(figsize=(W, 125 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1.0], hspace=0.7, wspace=0.75)
    ax, a1, a2 = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2])

    # A - total and minor-allele copy number along the genome
    for k, c in enumerate(chroms):
        if k % 2 == 0:
            ax.axvspan(offs[c], offs[c] + lens[c], color="#F3F3F3", lw=0, zorder=0)
    ymax = int(min(max(seg.CNt.max(), base + 2), 8))
    ax.axhline(base, ls=(0, (3, 2)), lw=0.7, color=GREY, zorder=1)
    for g0, g1, cn, b in seg[["g0", "g1", "CNt", "B"]].itertuples(index=False):
        ax.plot([g0, g1], [min(cn, ymax) + 0.08] * 2, color=INK, lw=2.2, solid_capstyle="butt", zorder=2)
        if cn > ymax:
            ax.plot((g0 + g1) / 2, ymax + 0.35, "^", ms=3, color=INK, clip_on=False)
        if b == b:
            ax.plot([g0, g1], [min(b, ymax) - 0.08] * 2, color=VERM, lw=1.8, solid_capstyle="butt", zorder=3)
    ax.set_xlim(0, lens.sum()); ax.set_ylim(-0.4, ymax + 0.5); ax.set_yticks(range(ymax + 1))
    ax.set_xticks((offs + lens / 2).values)
    ax.set_xticklabels([("\n" if c in ("20", "22") else "") + c for c in chroms], fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Copy number")
    ax.legend(handles=[Line2D([], [], color=INK, lw=2.2, label="Total copy number"),
                       Line2D([], [], color=VERM, lw=2.2, label="Minor allele (0 = LOH)"),
                       Line2D([], [], color=GREY, ls=(0, (3, 2)), lw=0.8, label=f"{base} copies (ploidy {ploidy:.2f})")],
              ncol=3, loc="lower right", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Allele-specific copy number", loc="left"); letter(ax, "A", x=-0.06)

    # B - arm-level gain / loss / LOH
    arms = arm_table(seg, base)
    x = np.arange(len(arms))
    a1.bar(x, arms.gain, color=VERM, width=0.75, edgecolor="none", label=f"Gain (> {base} copies)")
    a1.bar(x, -arms.loss, color=BLUE, width=0.75, edgecolor="none", label=f"Loss (< {base} copies)")
    loh = (arms.loh >= 50).values
    a1.scatter(x[loh], [110] * int(loh.sum()), marker="v", s=9, color=INK, zorder=3, label="LOH >= 50% of arm")
    a1.axhline(0, color=INK, lw=0.5)
    a1.set_xticks(x); a1.set_xticklabels(arms.arm, rotation=90, fontsize=5); a1.tick_params(axis="x", length=0)
    a1.set_xlim(-0.7, len(arms) - 0.3); a1.set_ylim(-108, 118)
    a1.set_yticks([-100, -50, 0, 50, 100]); a1.set_yticklabels(["100", "50", "0", "50", "100"])
    a1.set_ylabel("Arm affected (%)\n<- loss    gain ->")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, columnspacing=1.5)
    a1.set_title("Arm-level changes (autosomes)", loc="left"); letter(a1, "B", x=-0.09)

    # C - Sequenza's cellularity/ploidy solutions
    sc = a2.scatter(sol.ploidy, sol.cellularity, c=sol.SLPP, cmap="viridis", s=18, edgecolor="none", zorder=2)
    a2.scatter(sol.ploidy[:1], sol.cellularity[:1], s=80, facecolor="none", edgecolor=INK, lw=1.0, zorder=3)
    a2.set_xlabel("Ploidy"); a2.set_ylabel("Cellularity"); a2.set_ylim(0, 1.05)
    cb = fig.colorbar(sc, ax=a2, fraction=0.08, pad=0.04)
    cb.set_label("Scaled log posterior (SLPP)", fontsize=6.5); cb.ax.tick_params(labelsize=6, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)
    a2.set_title(f"Purity-ploidy solutions\nselected (circled): cellularity {cell:.2f}, ploidy {ploidy:.2f}", loc="left", fontsize=6.8); letter(a2, "C", x=-0.32)
    save(fig, a.o or "Figure4")

    print(f"\nSequenza selected solution: cellularity {cell:.2f}, ploidy {ploidy:.2f}; "
          f"gain/loss called against {base} copies")
    if len(sol) > 1:
        print("next solutions: " + ", ".join(f"{r.cellularity:.2f}/{r.ploidy:.2f} (SLPP {r.SLPP:.1f})"
                                              for r in sol.iloc[1:4].itertuples()))
    print(f"arms with LOH >= 50%: {', '.join(arms.arm[loh]) or 'none'}")
    ci = find(a.results, f"{a.tumor}_confints_CP.txt", required=False)
    if ci:
        print("Sequenza confidence-interval file (its point estimate should equal the selected solution):")
        print(open(ci).read().rstrip())
    print("Files used:\n  " + segf + "\n  " + solf + (f"\n  {ci}" if ci else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("qc", "cnv"):
        p = sp.add_parser(name)
        p.add_argument("--results", required=True, help="PAPIO somatic output folder")
        p.add_argument("--tumor", required=True, help="tumor sample ID, e.g. valid1_T")
        if name == "qc":
            p.add_argument("--normal", required=True, help="normal sample ID, e.g. valid1_N")
            p.add_argument("--min-relatedness", type=float, default=0.90,
                           help="tumor-normal relatedness threshold (pipeline: somalier_min_relatedness, 0.90)")
            p.add_argument("--max-ibs0", type=float, default=0.02,
                           help="IBS0/n at or below which low relatedness is LOH, not a swap (pipeline: somalier_max_ibs0_frac, 0.02)")
            p.add_argument("--max-contamination", type=float, default=0.03,
                           help="contamination threshold as a fraction (pipeline: qc_max_contamination, 0.03)")
        p.add_argument("-o", help="output name (default Figure2 / Figure4)")
    a = ap.parse_args()
    a.results = os.path.abspath(a.results)
    {"qc": qc, "cnv": cnv}[a.cmd](a)


if __name__ == "__main__":
    main()
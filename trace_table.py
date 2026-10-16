#!/usr/bin/env python3
"""
trace_table.py - per-process table of requested vs used CPU, memory and time
from a Nextflow trace file (Table 2 of the PAPIO protocol).

usage:
  python3 trace_table.py pipeline_trace.txt [--match TEXT] [-o Table2.tsv]
  --match keeps tasks whose name contains TEXT (e.g. a sample tag such as SRR7890883); omit it to use all tasks.

Per process (tasks that finished: COMPLETED or CACHED):
  CPU    : requested | used = %cpu/100 (median, max) | used/requested per task (median)
  Memory : requested (max; includes retries) | used max | used/requested per task (median, max)
  Time   : requested (needs the 'time' trace field) | wall time median, max
  CPU-hours used (wall x CPUs used) and reserved (wall x CPUs requested; what the cluster blocks)
  Suggested CPUs = busiest task rounded up; suggested memory = largest peak + 20%
  failed_attempts = FAILED tasks (later retried or not)

Memory used is peak_rss. If the trace has only rss (a sampled value, NOT the peak),
memory-used columns are labelled _rss and no memory is suggested, because rss can be
far below what a task needed. For a complete table add to nextflow.config:
  trace.fields = '...,peak_rss,peak_vmem,time,attempt'
"""
import argparse, math, re, sys
import pandas as pd


def die(msg):
    sys.exit(f"ERROR: {msg}")


def to_gb(s):
    """'16 GB', '1.2 GB', '500 MB', '1024 KB', '-' -> GB (Nextflow units, 2^30 bytes); NaN if absent"""
    m = re.match(r"\s*([\d.]+)\s*([KMGTP]?B)\s*$", str(s))
    if not m:
        return float("nan")
    return float(m.group(1)) * {"B": 2**-30, "KB": 2**-20, "MB": 2**-10, "GB": 1.0, "TB": 2**10, "PB": 2**20}[m.group(2)]


def to_min(s):
    """'1h 2m 3s', '45.2s', '120ms', '2d 3h' -> minutes; NaN if absent"""
    parts = re.findall(r"([\d.]+)\s*(ms|d|h|m|s)\b", str(s))
    if not parts:
        return float("nan")
    f = {"ms": 1 / 60000, "s": 1 / 60, "m": 1, "h": 60, "d": 1440}
    return sum(float(v) * f[u] for v, u in parts)


def to_pct(s):
    try:
        return float(str(s).rstrip("%"))
    except ValueError:
        return float("nan")


def ok(x):
    return x is not None and x == x


def pct(x):
    """ratio -> whole percent, half rounded up; '-' if missing"""
    return int(100 * x + 0.5) if ok(x) else "-"


def r1(x):
    return round(x, 1) if ok(x) else "-"


def fmt_time(minutes):
    if not ok(minutes):
        return "-"
    if minutes >= 90:
        return f"{minutes / 60:.1f} h"
    return f"{int(minutes + 0.5)} min" if minutes >= 1 else f"{int(minutes * 60 + 0.5)} s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace")
    ap.add_argument("--match", help="keep tasks whose name contains this text")
    ap.add_argument("-o", default="Table2.tsv")
    a = ap.parse_args()

    try:
        tr = pd.read_csv(a.trace, sep="\t", dtype=str)
    except FileNotFoundError:
        die(f"file not found: {a.trace}")
    need = {"name", "status", "realtime", "cpus", "memory", "%cpu"}
    if need - set(tr.columns):
        die(f"trace lacks columns: {', '.join(sorted(need - set(tr.columns)))}")
    if a.match:
        tr = tr[tr["name"].str.contains(a.match, regex=False)]
        if tr.empty:
            die(f"no task name contains '{a.match}'; run without --match or use a tag from the name column")

    peak = "peak_rss" in tr.columns
    mem_col = "peak_rss" if peak else "rss"
    if mem_col not in tr.columns:
        die("trace has neither peak_rss nor rss")
    has_time = "time" in tr.columns
    notes = []
    if not peak:
        notes.append("no peak_rss in the trace: memory used = rss (sampled, can be far below the peak); "
                     "no memory suggestion made")
    if not has_time:
        notes.append("no 'time' field in the trace: requested time not shown")

    base = tr["name"].str.replace(r"\s*\(.*\)\s*$", "", regex=True)
    tr["process"] = base.str.split(":").str[-1]
    paths = base.groupby(tr["process"]).nunique()
    for p in paths[paths > 1].index:
        notes.append(f"{p}: {paths[p]} workflow paths share this process name and are merged in one row")

    tr["cpu_req"] = pd.to_numeric(tr["cpus"], errors="coerce")
    tr["cpu_used"] = tr["%cpu"].map(to_pct) / 100
    tr["mem_req"] = tr["memory"].map(to_gb)
    tr["mem_used"] = tr[mem_col].map(to_gb)
    tr["wall_min"] = tr["realtime"].map(to_min)
    tr["time_req"] = tr["time"].map(to_min) if has_time else float("nan")
    tr["cpu_ratio"] = tr["cpu_used"] / tr["cpu_req"]          # per task, so retries with bigger
    tr["mem_ratio"] = tr["mem_used"] / tr["mem_req"]          # requests are compared like for like

    done = tr[tr["status"].isin(["COMPLETED", "CACHED"])]
    failed = tr[tr["status"] == "FAILED"].groupby("process").size()
    if done.empty:
        die("no COMPLETED or CACHED tasks in the trace")

    mu = "mem_used_max_GB" if peak else "mem_rss_max_GB"
    rows = []
    for p in sorted(set(tr["process"])):                     # includes processes that only failed
        g = done[done["process"] == p]
        cmax, mmax = g["cpu_used"].max(), g["mem_used"].max()
        rows.append({
            "process": p,
            "tasks": len(g),
            "failed_attempts": int(failed.get(p, 0)),
            "cpus_requested": int(g["cpu_req"].max()) if ok(g["cpu_req"].max()) else "-",
            "cpus_used_median": r1(g["cpu_used"].median()),
            "cpus_used_max": r1(cmax),
            "cpu_use_%": pct(g["cpu_ratio"].median()),
            "mem_requested_GB": r1(g["mem_req"].max()),
            mu: r1(mmax),
            "mem_use_%_median": pct(g["mem_ratio"].median()),
            "mem_use_%_max": pct(g["mem_ratio"].max()),
            "time_requested": fmt_time(g["time_req"].max()),
            "wall_median": fmt_time(g["wall_min"].median()),
            "wall_max": fmt_time(g["wall_min"].max()),
            "cpu_hours_used": round((g["wall_min"] / 60 * g["cpu_used"]).sum(), 1),
            "cpu_hours_reserved": round((g["wall_min"] / 60 * g["cpu_req"]).sum(), 1),
            "suggested_cpus": max(1, math.ceil(cmax - 0.05)) if ok(cmax) else "-",
            "suggested_mem_GB": math.ceil(mmax * 1.2) if peak and ok(mmax) else "-",
        })
    out = pd.DataFrame(rows).sort_values("cpu_hours_reserved", ascending=False, kind="stable")
    if not has_time:
        out = out.drop(columns="time_requested")

    total = {c: "" for c in out.columns}
    total.update(process="TOTAL", tasks=int(out["tasks"].sum()), failed_attempts=int(out["failed_attempts"].sum()),
                 cpu_hours_used=round((done["wall_min"] / 60 * done["cpu_used"]).sum(), 1),
                 cpu_hours_reserved=round((done["wall_min"] / 60 * done["cpu_req"]).sum(), 1))
    out = pd.concat([out, pd.DataFrame([total])], ignore_index=True)
    out.to_csv(a.o, sep="\t", index=False)

    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(out.to_string(index=False))
    if "submit" in tr.columns and "complete" in tr.columns:
        s = pd.to_datetime(tr["submit"], errors="coerce").min()
        e = pd.to_datetime(tr["complete"], errors="coerce").max()
        if ok(s) and ok(e):
            print(f"\nfirst submit to last completion: {(e - s).total_seconds() / 3600:.1f} h "
                  "(includes queueing; with -resume, CACHED tasks keep their original times)")
    for n in notes:
        print("note:", n)
    print(f"wrote {a.o}")


if __name__ == "__main__":
    main()
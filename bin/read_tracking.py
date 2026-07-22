#!/usr/bin/env python3
"""Per-sample read accounting, starting from the raw FASTQ.

Every downstream stats file starts *after* primer removal, so a sample whose
reads were discarded by cutadapt (wrong primer pair for that run) looked simply
"small" rather than "trimmed to death" — the raw count is the only place the
distinction shows up. Aggregate the cutadapt logs (and the DADA2 filter stats
when present) into one table with a `reads_raw` column plus a TOTAL row.

Usage: read_tracking.py <out.tsv> [cutadapt_logs...] [filt_stats...]
"""
import glob
import os
import re
import sys

RAW_RE = re.compile(r"Total read pairs processed:\s*([\d,]+)")
RAW_SE_RE = re.compile(r"Total reads processed:\s*([\d,]+)")
KEPT_RE = re.compile(r"Pairs written \(passing filters\):\s*([\d,]+)")
KEPT_SE_RE = re.compile(r"Reads written \(passing filters\):\s*([\d,]+)")


def _num(m):
    return int(m.group(1).replace(",", "")) if m else 0


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: read_tracking.py <out.tsv> [inputs...]")
    out_path = sys.argv[1]
    inputs = sys.argv[2:] or sorted(glob.glob("*_cutadapt.log")) + sorted(glob.glob("*_filt_stats.tsv"))

    samples = {}

    for f in inputs:
        base = os.path.basename(f)
        if base.endswith("_cutadapt.log"):
            acc = base[: -len("_cutadapt.log")]
            try:
                txt = open(f, errors="ignore").read()
            except OSError:
                continue
            raw = _num(RAW_RE.search(txt)) or _num(RAW_SE_RE.search(txt))
            kept = _num(KEPT_RE.search(txt)) or _num(KEPT_SE_RE.search(txt))
            d = samples.setdefault(acc, {})
            d["reads_raw"] = raw
            d["reads_after_primer"] = kept
        elif base.endswith("_filt_stats.tsv"):
            acc = base[: -len("_filt_stats.tsv")]
            try:
                lines = [l.rstrip("\n") for l in open(f) if l.strip()]
            except OSError:
                continue
            if len(lines) < 2:
                continue
            row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            d = samples.setdefault(acc, {})
            try:
                d["reads_after_filter"] = int(float(row.get("reads_out", 0)))
            except ValueError:
                pass

    cols = ["reads_raw", "reads_after_primer", "reads_after_filter"]
    with open(out_path, "w") as out:
        out.write("sample\t" + "\t".join(cols) + "\tpct_primer\tpct_overall\n")
        totals = dict.fromkeys(cols, 0)
        for acc in sorted(samples):
            d = samples[acc]
            vals = [d.get(c, 0) for c in cols]
            for c, v in zip(cols, vals):
                totals[c] += v
            raw = vals[0] or 1
            out.write(
                f"{acc}\t" + "\t".join(str(v) for v in vals)
                + f"\t{100 * vals[1] / raw:.1f}\t{100 * vals[2] / raw:.1f}\n"
            )
        raw = totals["reads_raw"] or 1
        out.write(
            "TOTAL\t" + "\t".join(str(totals[c]) for c in cols)
            + f"\t{100 * totals['reads_after_primer'] / raw:.1f}"
            + f"\t{100 * totals['reads_after_filter'] / raw:.1f}\n"
        )
    # Also emit JSON next to the TSV for the viz Provenance tab. Written as
    # data/provenance.json when placed in the viz payload.
    import json
    json_path = os.path.join(os.path.dirname(out_path) or ".", "provenance.json")
    stages = [
        ("raw", "Raw FASTQ", "reads_raw"),
        ("primer", "After primer removal", "reads_after_primer"),
        ("filter", "After quality filter", "reads_after_filter"),
    ]
    payload = {
        "stages": [{"id": sid, "label": lab} for sid, lab, _ in stages],
        "total": {sid: totals[col] for sid, _, col in stages},
        "samples": {
            acc: {sid: samples[acc].get(col, 0) for sid, _, col in stages}
            for acc in sorted(samples)
        },
    }
    try:
        with open(json_path, "w") as jf:
            json.dump(payload, jf)
    except OSError as e:
        print(f"[WARN] could not write {json_path}: {e}")

    print(f"[INFO] read_tracking: {len(samples)} samples, "
          f"{totals['reads_raw']:,} raw reads -> "
          f"{totals['reads_after_primer']:,} after primers")


if __name__ == "__main__":
    main()

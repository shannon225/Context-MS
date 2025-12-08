#!/usr/bin/env python3
"""
Count peptides per seed passing pyIsoPEP thresholds (PEP and q),
deduplicating by 'id' (keep highest 'primary'), for multiple thresholds.
Defaults to thresholds 0.01 and 0.05.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd

def pick_col(df, preferred_exact, predicate):
    for name in preferred_exact:
        if name in df.columns:
            return name
    for name in preferred_exact:
        for c in df.columns:
            if c.lower() == name.lower():
                return c
    matches = [c for c in df.columns if predicate(c)]
    return matches[0] if matches else None

def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    df.columns = df.columns.str.strip()
    return df

def extract_seed(path: Path) -> str:
    m = re.search(r"seed(\d+)", path.name)
    return m.group(1) if m else path.stem

def deduplicate_by_id(df: pd.DataFrame) -> pd.DataFrame:
    if "id" not in df.columns or "primary" not in df.columns:
        return df
    df["primary"] = pd.to_numeric(df["primary"], errors="coerce").fillna(float("-inf"))
    df = (df.sort_values(["id", "primary"], ascending=[True, False])
            .drop_duplicates("id", keep="first"))
    return df

def main():
    ap = argparse.ArgumentParser(
        description="Count pyIsoPEP PEP/q passing thresholds per seed (dedup by id)."
    )
    ap.add_argument("--indir", type=Path, default=Path("3_pep"),
                    help="Directory containing *_pep_seed*_ipspline.txt")
    ap.add_argument("--pattern", default="*_pep_seed*_ipspline.txt",
                    help="Glob pattern inside --indir to find PEP tables")
    ap.add_argument("--thresh", type=float, nargs="+",
                    help="One or more thresholds (default: 0.01 and 0.05)")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output TSV path (default: 3_pep/pyispep_counts_multi_dedup.tsv)")
    args = ap.parse_args()

    # Default thresholds if not provided
    if not args.thresh:
        args.thresh = [0.01, 0.05]

    files = sorted(Path(args.indir).glob(args.pattern))
    if not files:
        sys.exit(f"No files matched: {args.indir}/{args.pattern}")

    rows = []
    for f in files:
        try:
            df = load_table(f)
        except Exception as e:
            print(f"WARN: failed to read {f}: {e}", file=sys.stderr)
            continue

        before = len(df)
        df = deduplicate_by_id(df)
        after = len(df)
        if after < before:
            print(f"[{f.name}] Deduplicated {before} → {after} rows", file=sys.stderr)

        pep_col = pick_col(df, ["pyIsoPEP PEP"],
                           lambda c: ("pyispep" in c.lower() and "pep" in c.lower() and "q" not in c.lower()))
        q_col   = pick_col(df, ["pyIsoPEP q-value from FDR"],
                           lambda c: ("pyispep" in c.lower() and "q" in c.lower() and "fdr" in c.lower()))
        label_c = next((c for c in df.columns if c.lower() == "label"), None)

        if not pep_col or not q_col:
            print(f"WARN: missing PEP/q columns in {f.name} (pep_col={pep_col}, q_col={q_col})", file=sys.stderr)
            continue

        df[pep_col] = pd.to_numeric(df[pep_col], errors="coerce")
        df[q_col]   = pd.to_numeric(df[q_col], errors="coerce")

        targ = None
        if label_c:
            lab = pd.to_numeric(df[label_c].astype(str).str.strip(), errors="coerce").fillna(0).astype(int)
            targ = (lab == 1)

        total_rows = int(len(df))
        seed = extract_seed(f)

        for thr in args.thresh:
            pep_ok = df[pep_col] < thr
            q_ok   = df[q_col]   < thr
            both   = pep_ok & q_ok

            rec = {
                "seed": seed,
                "file": f.name,
                "threshold": thr,
                "total_rows": total_rows,
                "pep_lt": int(pep_ok.sum()),
                "qval_lt": int(q_ok.sum()),
                "both_lt": int(both.sum()),
            }
            if targ is not None:
                rec.update({
                    "targets_rows": int(targ.sum()),
                    "targets_pep_lt": int((targ & pep_ok).sum()),
                    "targets_qval_lt": int((targ & q_ok).sum()),
                    "targets_both_lt": int((targ & both).sum()),
                })
            rows.append(rec)

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        sys.exit("No valid data to summarize.")

    out_df["seed_num"] = pd.to_numeric(out_df["seed"], errors="coerce")
    out_df = (out_df.sort_values(["seed_num", "seed", "threshold"], na_position="last")
                      .drop(columns=["seed_num"]))

    out_path = args.out or (args.indir / "pyispep_counts_multi_dedup.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"✅ Wrote summary: {out_path}")
    print(out_df.head(min(12, len(out_df))).to_string(index=False))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Aggregate per-chromosome cell counts from __1_PerChrVariantJoin.py into the
gnomAD-vs-LRS variant-density enrichment table.

Enrichment definition mirrors the gnomAD contingency analysis exactly
(matching __1_+__10_ContingencyStatistics):
  enr = (a / nPos_GC) / (c / nPos_nonGC)
where the partition is hierarchical/cumulative:
  GC at k=K     = positions with longest_k = K (exclusive bucket)
  non-GC at k=K = positions with longest_k < K + pure-nonGC
                  (longer-k positions are tested in their own k cell)

Position denominators (a+b, c+d) are taken directly from the
StatsSummary.R.out.txt files produced by step 06+07, ensuring exact identity
with the published gnomAD enrichment.

The same formula and denominators are applied with GA4K SNVs as variant
source (enr_lrs) — this is the orthogonal HiFi reproduction of the gnomAD
enrichment.

Output:
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_enrichment_compare.tsv

Environment:
  PROJECT_ROOT must point to the project root.
"""
import collections, glob, math, csv, os, sys
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", ""))
if not PROJECT_ROOT.is_dir():
    sys.exit("PROJECT_ROOT env var must point to a directory")

BASE    = PROJECT_ROOT / "geneconv_complete"
ART     = BASE / "artefact_validation_ga4k"
IN_DIR  = ART / "results" / "per_chr_v2"
OUT_DIR = ART / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Aggregate per-chr cell counts
agg = collections.defaultdict(lambda: [0, 0, 0.0, 0])  # n_ann, n_conf, sum_logr, n_lrs
nfiles = 0
for f in sorted(glob.glob(str(IN_DIR / "chr*.cell_counts.tsv"))):
    nfiles += 1
    with open(f) as h:
        header = next(h)
        if "n_lrs" not in header:
            sys.exit(f"ERROR: {f} missing n_lrs column - rerun worker first")
        for line in h:
            p = line.rstrip().split("\t")
            if len(p) < 7: continue
            s, k, ab = p[0], p[1], p[2]
            cell = agg[(s, k, ab)]
            cell[0] += int(p[3])
            cell[1] += int(p[4])
            cell[2] += float(p[5])
            cell[3] += int(p[6])
print(f"aggregated {len(agg)} cells from {nfiles} chrs")

# 2. Load ST3 position denominators (nPos_GC = a+b, nPos_nonGC = c+d per cell)
STRATA = ["allmapp", "segdupmapp", "nosegdupmapp"]
AF_LBL_TO_CODE = {
    "0.001_0.005": "0001_0005",
    "0.005_0.01":  "0005_001",
    "0.01_0.05":   "001_005",
    "0.05_0.1":    "005_01",
    "0.1_0.5":     "01_05",
    "0.5_2":       "05_2",
}
nPos = {}
for s in STRATA:
    fp = BASE / f"stats_{s}" / "_All_1000000_1000000_gnomad~genome_StatsSummary.R.out.txt"
    if not fp.exists():
        print(f"WARN missing stats file: {fp}"); continue
    with open(fp) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = row["RepLength"]
            af_code = row["FrequencyInterval"]
            af_lbl = next((lbl for lbl, c in AF_LBL_TO_CODE.items() if c == af_code), None)
            if af_lbl is None: continue
            try:
                a = int(row["NrOfVarPosConvPos"])
                b = int(row["NrOfNoVarPosConvPos"])
                c = int(row["NrOfVarPosNoConvPos"])
                d = int(row["NrOfNoVarPosNoConvPos"])
            except (KeyError, ValueError): continue
            nPos[(s, k, af_lbl)] = (a + b, c + d)
print(f"loaded {len(nPos)} position-denominator cells from ST3 stats")

# 3. Compute both enrichments per cell with hierarchical partition
KS = ["17","19","21","31","41","51","61","71","81","91"]
ALL_CLASSES = KS + ["nonGC"]
AFS = ["0.001_0.005","0.005_0.01","0.01_0.05","0.05_0.1","0.1_0.5","0.5_2"]

out_path = OUT_DIR / "v2_enrichment_compare.tsv"
with open(out_path, "w") as o:
    o.write("stratum\tk\taf_bin\t"
            "a_gnomad\tc_gnomad\tnPos_GC\tnPos_nonGC\t"
            "a_lrs\tc_lrs\t"
            "enr_gnomad\tenr_lrs\tratio_lrs_over_gnomad\n")
    print(f"\n{'stratum':<13} {'k':>3} {'af':>11} "
          f"{'enr_gnomad':>11} {'enr_lrs':>9} {'lrs/gn':>7}")
    for s in STRATA:
        for k in KS:
            for ab in AFS:
                gc = agg.get((s, k, ab), [0,0,0.0,0])
                # Hierarchical partition: non-GC = shorter-k + pure-nonGC
                a_lrs_ng = 0; a_gn_ng = 0
                for kc in ALL_CLASSES:
                    if kc == k: continue
                    if kc != "nonGC" and int(kc) > int(k): continue
                    cell = agg.get((s, kc, ab), [0,0,0.0,0])
                    a_gn_ng += cell[0]
                    a_lrs_ng += cell[3]
                ab_count, cd_count = nPos.get((s, k, ab), (0, 0))
                if ab_count == 0 or cd_count == 0: continue
                a_gn = gc[0]; c_gn = a_gn_ng
                a_lr = gc[3]; c_lr = a_lrs_ng
                dgc_gn = a_gn / ab_count
                dnc_gn = c_gn / cd_count
                enr_gn = (dgc_gn / dnc_gn) if dnc_gn > 0 else float("nan")
                dgc_lr = a_lr / ab_count
                dnc_lr = c_lr / cd_count
                enr_lr = (dgc_lr / dnc_lr) if dnc_lr > 0 else float("nan")
                ratio = (enr_lr / enr_gn) if (enr_gn and not math.isnan(enr_gn)) else float("nan")
                o.write(f"{s}\t{k}\t{ab}\t{a_gn}\t{c_gn}\t{ab_count}\t{cd_count}\t"
                        f"{a_lr}\t{c_lr}\t{enr_gn:.4f}\t{enr_lr:.4f}\t{ratio:.4f}\n")
                if s == "segdupmapp" and int(k) >= 51:
                    print(f"{s:<13} {k:>3} {ab:>11} "
                          f"{enr_gn:>11.3f} {enr_lr:>9.3f} {ratio:>7.3f}")
        if s == "segdupmapp": print()
print(f"\nwrote {out_path}")

#!/usr/bin/env python3
"""
Inverse-variance-pooled risk difference of HiFi mismatch rate
(1 - confirmation rate) between GC and non-GC variants per (mappability
stratum, k), pooled across the 6 AF strata, with hierarchical partition
matching __1_+__10_ContingencyStatistics.

Outputs per (stratum, k):
  rd_pooled    AF-adjusted Delta-mismatch (GC - non-GC)
  ci_lo, ci_hi 95% confidence interval
  max_inflation_at_ci_hi = 1 / (1 - ci_hi)  -- worst-case enrichment inflation
                                                attributable to mismapping at upper CI

Output:
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_mh_pooled.tsv
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_mh_per_af.tsv

Environment:
  PROJECT_ROOT must point to the project root.
"""
import collections, glob, math, os, sys
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", ""))
if not PROJECT_ROOT.is_dir():
    sys.exit("PROJECT_ROOT env var must point to a directory")

ART     = PROJECT_ROOT / "geneconv_complete" / "artefact_validation_ga4k"
IN_DIR  = ART / "results" / "per_chr_v2"
OUT_DIR = ART / "results"

agg = collections.defaultdict(lambda: [0, 0])  # n_annotated, n_confirmed
for f in sorted(glob.glob(str(IN_DIR / "chr*.cell_counts.tsv"))):
    with open(f) as h:
        next(h)
        for line in h:
            p = line.rstrip().split("\t")
            if len(p) < 6: continue
            s, k, ab, n_ann, n_conf = p[0], p[1], p[2], p[3], p[4]
            c = agg[(s, k, ab)]
            c[0] += int(n_ann); c[1] += int(n_conf)

STRATA = ["allmapp", "segdupmapp", "nosegdupmapp"]
KS     = ["17","19","21","31","41","51","61","71","81","91"]
ALL_CLASSES = KS + ["nonGC"]
af_bins = sorted({k[2] for k in agg.keys()})
print(f"discovered {len(af_bins)} AF bins")

Z = 1.96
def pooled_rd(per_strat):
    """Inverse-variance pooled risk difference + 95% CI."""
    sum_w_rd = 0.0; sum_w = 0.0; n_used = 0; total_n = 0
    for (n1, a1, n0, a0) in per_strat:
        if n1 < 5 or n0 < 5: continue
        p1, p0 = a1/n1, a0/n0
        var = p1*(1-p1)/n1 + p0*(1-p0)/n0
        if var <= 0: continue
        w = 1/var
        sum_w_rd += w * (p1 - p0)
        sum_w   += w
        n_used  += 1
        total_n += n1 + n0
    if sum_w == 0: return (None, None, None, 0, 0)
    rd = sum_w_rd / sum_w
    se = math.sqrt(1/sum_w)
    return (rd, rd - Z*se, rd + Z*se, n_used, total_n)

# Per-AF detail table
out_per_af = OUT_DIR / "v2_mh_per_af.tsv"
with open(out_per_af, "w") as o:
    o.write("stratum\tk\taf_bin\tn_gc\tmis_gc\tn_nongc\tmis_nongc\trd\tci_lo\tci_hi\n")
    for s in STRATA:
        for k in KS:
            for ab in af_bins:
                gc = agg.get((s, k, ab), [0,0])
                # Hierarchical partition
                ng_a = ng_c = 0
                for kc in ALL_CLASSES:
                    if kc == k: continue
                    if kc != "nonGC" and int(kc) > int(k): continue
                    cell = agg.get((s, kc, ab), [0,0])
                    ng_a += cell[0]; ng_c += cell[1]
                gc_a, gc_c = gc[0], gc[1]
                if gc_a < 5 or ng_a < 5: continue
                gc_u = gc_a - gc_c; ng_u = ng_a - ng_c
                p1 = gc_u/gc_a; p0 = ng_u/ng_a
                rd = p1 - p0
                var = p1*(1-p1)/gc_a + p0*(1-p0)/ng_a
                se = math.sqrt(var) if var > 0 else 0
                lo, hi = rd - Z*se, rd + Z*se
                o.write(f"{s}\t{k}\t{ab}\t{gc_a}\t{p1:.4f}\t"
                        f"{ng_a}\t{p0:.4f}\t{rd:+.4f}\t{lo:+.4f}\t{hi:+.4f}\n")
print(f"wrote {out_per_af}")

# Pooled-per-(stratum, k) table + max enrichment-inflation bound
out_pooled = OUT_DIR / "v2_mh_pooled.tsv"
with open(out_pooled, "w") as o:
    o.write("stratum\tk\tn_strata_used\ttotal_n\trd_pooled\tci_lo\tci_hi\t"
            "max_inflation_at_ci_hi\n")
    print(f"\n{'stratum':<13} {'k':>3} {'n_strata':>9} {'total_n':>12} "
          f"{'RD_pooled':>10} {'95% CI':>22} {'max_infl':>10}")
    for s in STRATA:
        for k in KS:
            per_strat = []
            for ab in af_bins:
                gc = agg.get((s, k, ab), [0,0])
                ng_a = ng_c = 0
                for kc in ALL_CLASSES:
                    if kc == k: continue
                    if kc != "nonGC" and int(kc) > int(k): continue
                    cell = agg.get((s, kc, ab), [0,0])
                    ng_a += cell[0]; ng_c += cell[1]
                gc_a, gc_c = gc[0], gc[1]
                if gc_a < 5 or ng_a < 5: continue
                per_strat.append((gc_a, gc_a - gc_c, ng_a, ng_a - ng_c))
            rd, lo, hi, n_used, ttl = pooled_rd(per_strat)
            if rd is None: continue
            max_infl = 1/(1-hi) if hi < 1 else float("inf")
            ci_str = f"[{lo:+.4f}, {hi:+.4f}]"
            o.write(f"{s}\t{k}\t{n_used}\t{ttl}\t{rd:+.4f}\t{lo:+.4f}\t{hi:+.4f}\t{max_infl:.3f}\n")
            print(f"{s:<13} {k:>3} {n_used:>9} {ttl:>12} "
                  f"{rd:>+10.4f} {ci_str:>22} {max_infl:>10.3f}")
print(f"\nwrote {out_pooled}")

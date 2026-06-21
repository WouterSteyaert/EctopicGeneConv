#!/usr/bin/env python3
"""
Pooled HiFi confirmation rates for GC and non-GC variants per (mappability
stratum, k), restricted to the cohort-saturated reliable AF zone (AF >= 0.01).
Hierarchical partition matches __1_+__10_ContingencyStatistics.

Wilson 95% confidence interval on each rate.

Output (used by __3_MakeConfirmationFigure.R for the headline figure):
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_conf_rates_reliableAF.tsv

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

agg = collections.defaultdict(lambda: [0, 0, 0.0, 0])
for f in sorted(glob.glob(str(IN_DIR / "chr*.cell_counts.tsv"))):
    with open(f) as h:
        next(h)
        for line in h:
            p = line.rstrip().split("\t")
            if len(p) < 6: continue
            s, k, ab = p[0], p[1], p[2]
            cell = agg[(s, k, ab)]
            cell[0] += int(p[3])  # n_annotated
            cell[1] += int(p[4])  # n_confirmed

STRATA = ["allmapp", "segdupmapp", "nosegdupmapp"]
KS     = ["17","19","21","31","41","51","61","71","81","91"]
ALL_CLASSES = KS + ["nonGC"]
RELIABLE = ["0.01_0.05","0.05_0.1","0.1_0.5","0.5_2"]   # AF >= 0.01

Z = 1.96
def wilson(k, n):
    if n == 0: return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + Z*Z/n
    c = (p + Z*Z/(2*n)) / den
    h = Z * math.sqrt(p*(1-p)/n + Z*Z/(4*n*n)) / den
    return (p, max(0, c-h), min(1, c+h))

out = OUT_DIR / "v2_conf_rates_reliableAF.tsv"
with open(out, "w") as o:
    o.write("stratum\tk\tn_GC_ann\tconf_GC\tGC_lo\tGC_hi\t"
            "n_nonGC_ann\tconf_nonGC\tnonGC_lo\tnonGC_hi\n")
    print(f"{'stratum':<13} {'k':>3} "
          f"{'n_GC':>9} {'GC':>17}  "
          f"{'n_nonGC':>11} {'nonGC':>17}")
    for s in STRATA:
        for k in KS:
            gc_a = gc_c = 0
            ng_a = ng_c = 0
            for ab in RELIABLE:
                gc_cell = agg.get((s, k, ab), [0,0,0,0])
                gc_a += gc_cell[0]; gc_c += gc_cell[1]
                # Hierarchical partition
                for kc in ALL_CLASSES:
                    if kc == k: continue
                    if kc != "nonGC" and int(kc) > int(k): continue
                    ng_cell = agg.get((s, kc, ab), [0,0,0,0])
                    ng_a += ng_cell[0]; ng_c += ng_cell[1]
            if gc_a < 5 or ng_a < 5: continue
            mg, mg_lo, mg_hi = wilson(gc_c, gc_a)
            mn, mn_lo, mn_hi = wilson(ng_c, ng_a)
            o.write(f"{s}\t{k}\t{gc_a}\t{mg:.4f}\t{mg_lo:.4f}\t{mg_hi:.4f}\t"
                    f"{ng_a}\t{mn:.4f}\t{mn_lo:.4f}\t{mn_hi:.4f}\n")
            gc_str = f"{mg:.4f}[{mg_lo:.4f},{mg_hi:.4f}]"
            ng_str = f"{mn:.4f}[{mn_lo:.4f},{mn_hi:.4f}]"
            print(f"{s:<13} {k:>3} {gc_a:>9} {gc_str:>17}  "
                  f"{ng_a:>11} {ng_str:>17}")
        print()
print(f"\nwrote {out}")

#!/usr/bin/env python3
"""
Per-cell reference aggregation of __1_PerChrVariantJoin.py outputs.

Produces the per-cell raw counts and per-AF GC vs non-GC contrast tables
(reference data; the headline tables come from __2_AggregateMhRd.py,
__2_AggregateConfRates.py, and __2_AggregateEnrichment.py).

Hierarchical partition matches __1_+__10_ContingencyStatistics.

Outputs:
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_per_cell_aggregated.tsv
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_contrast_per_af.tsv

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

# Per-cell: n_ann, n_conf, sum_log_ratio, n_logr (allele pair count)
agg = collections.defaultdict(lambda: [0, 0, 0.0, 0])
for f in sorted(glob.glob(str(IN_DIR / "chr*.cell_counts.tsv"))):
    with open(f) as h:
        next(h)
        for line in h:
            p = line.rstrip().split("\t")
            if len(p) < 6: continue
            s, k, ab, n_ann, n_conf, slr = p[:6]
            n_logr_v = int(p[7]) if len(p) >= 8 else int(n_conf)
            cell = agg[(s, k, ab)]
            cell[0] += int(n_ann); cell[1] += int(n_conf); cell[2] += float(slr)
            cell[3] += n_logr_v
print(f"aggregated {len(agg)} cells from {len(glob.glob(str(IN_DIR / 'chr*.cell_counts.tsv')))} chrs")

out_cells = OUT_DIR / "v2_per_cell_aggregated.tsv"
with open(out_cells, "w") as o:
    o.write("stratum\tk_class\taf_bin\tn_annotated\tn_confirmed\tconf_rate\tmean_log_ratio\n")
    for key in sorted(agg.keys()):
        s, k, ab = key
        n_a, n_c, slr, n_lr = agg[key]
        cr  = n_c / n_a if n_a else 0
        mlr = slr / n_lr if n_lr else 0
        o.write(f"{s}\t{k}\t{ab}\t{n_a}\t{n_c}\t{cr:.4f}\t{mlr:+.4f}\n")
print(f"wrote {out_cells}")

STRATA = ["allmapp", "segdupmapp", "nosegdupmapp"]
KS     = ["17","19","21","31","41","51","61","71","81","91"]
ALL_CLASSES = KS + ["nonGC"]
AFS    = ["0.001_0.005","0.005_0.01","0.01_0.05","0.05_0.1","0.1_0.5","0.5_2"]

out_c = OUT_DIR / "v2_contrast_per_af.tsv"
with open(out_c, "w") as o:
    o.write("stratum\tk\taf_bin\tn_gc_ann\tn_gc_conf\tconf_rate_gc\t"
            "n_nongc_ann\tn_nongc_conf\tconf_rate_nongc\t"
            "delta_conf_rate\tmean_lr_gc\tmean_lr_nongc\tdelta_log_ratio\n")
    print(f"\n{'stratum':<13} {'k':>3} {'AF':>11} "
          f"{'cr_gc':>7} {'cr_ng':>7} {'D_cr':>7} "
          f"{'lr_gc':>7} {'lr_ng':>7} {'D_lr':>7}")
    for s in STRATA:
        for k in KS:
            for ab in AFS:
                gc = agg.get((s, k, ab), [0,0,0.0,0])
                # Hierarchical partition: shorter-k + pure-nonGC
                ng_a = ng_c = 0; ng_slr = 0.0; ng_logr = 0
                for kc in ALL_CLASSES:
                    if kc == k: continue
                    if kc != "nonGC" and int(kc) > int(k): continue
                    cell = agg.get((s, kc, ab), [0,0,0.0,0])
                    ng_a += cell[0]; ng_c += cell[1]; ng_slr += cell[2]; ng_logr += cell[3]
                if not gc[0] or not ng_a: continue
                cr_g = gc[1]/gc[0]
                cr_n = ng_c/ng_a
                lr_g = gc[2]/gc[3] if gc[3] else 0
                lr_n = ng_slr/ng_logr if ng_logr else 0
                o.write(f"{s}\t{k}\t{ab}\t{gc[0]}\t{gc[1]}\t{cr_g:.4f}\t"
                        f"{ng_a}\t{ng_c}\t{cr_n:.4f}\t"
                        f"{cr_g-cr_n:+.4f}\t{lr_g:+.4f}\t{lr_n:+.4f}\t{lr_g-lr_n:+.4f}\n")
                print(f"{s:<13} {k:>3} {ab:>11} "
                      f"{cr_g:>7.3f} {cr_n:>7.3f} {cr_g-cr_n:>+7.3f} "
                      f"{lr_g:>+7.3f} {lr_n:>+7.3f} {lr_g-lr_n:>+7.3f}")
print(f"\nwrote {out_c}")

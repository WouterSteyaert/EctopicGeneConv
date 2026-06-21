#!/usr/bin/env python3
"""
Per-chromosome variant-level join: gnomAD genome SNVs x GA4K PacBio HiFi
joint VCF, classified by the set-difference longest-template-length k assignment
(matching __1_+__10_ContingencyStatistics partition exactly), mappability
stratum, and AF bin.

Counting is per POSITION (matches the perl `$Variants{$VarPos}` logic), so a
multi-allelic position with two alleles in the same AF bin counts once for
that bin (and once for each AF bin it spans).

Per (mappability_stratum, k_class, AF_bin) cell, the worker outputs:
  n_annotated   gnomAD positions in the cell
  n_confirmed   gnomAD positions where >= 1 (REF, ALT) matches a GA4K SNV
  sum_log_ratio Sigma log10(gnomAD_AF / GA4K_AF) summed over confirmed allele pairs
  n_lrs         GA4K HiFi positions in the cell (regardless of gnomAD)
  n_logr        count of allele pairs contributing to sum_log_ratio

Output (one per chromosome):
  <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/per_chr_v2/chr<CHR>.cell_counts.tsv

Usage:
  python3 __1_PerChrVariantJoin.py --chr=<1..22|X|Y>

Environment:
  PROJECT_ROOT must point to the project root (see 00_Configuration/README.md).
"""
import gzip, re, os, sys, collections, bisect, math, argparse
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--chr", required=True, help="Chromosome label (1..22, X, Y; no 'chr' prefix)")
args = ap.parse_args()
CHR = args.chr

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", ""))
if not PROJECT_ROOT.is_dir():
    sys.exit("PROJECT_ROOT env var must point to a directory")

BASE     = PROJECT_ROOT / "geneconv_complete"
GCDIFF   = BASE / "gcdiffs"
CONV91   = BASE / "conv" / "91"
REGIONS  = BASE / "regions"
GNOMAD   = BASE / "gnomad" / "genome" / f"{CHR}.sorted.bed.gz"
ART_BASE = BASE / "artefact_validation_ga4k"
GA4K_VCF = ART_BASE / "vcf" / f"pb_joint_merged.snv.chr{CHR}.vcf.gz"
OUT_DIR  = ART_BASE / "results" / "per_chr_v2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Position -> longest_k assignment (set-difference, matches __1_ perl)
K_SOURCES = [
    (17, GCDIFF / "17_19" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (19, GCDIFF / "19_21" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (21, GCDIFF / "21_31" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (31, GCDIFF / "31_41" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (41, GCDIFF / "41_51" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (51, GCDIFF / "51_61" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (61, GCDIFF / "61_71" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (71, GCDIFF / "71_81" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (81, GCDIFF / "81_91" / f"{CHR}.positions.gcdiffs.sorted.bed.gz"),
    (91, CONV91 / f"{CHR}.repsum.geneconv.sorted.bed.gz"),
]
pos_to_k = {}
for k, fp in K_SOURCES:
    if not fp.exists(): continue
    with gzip.open(fp, "rt") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) < 3: continue
            try: pos = int(cols[2])
            except: continue
            if pos not in pos_to_k:
                pos_to_k[pos] = k
print(f"chr{CHR}: {len(pos_to_k)} conv positions loaded", flush=True)

# 2. Mappability stratum interval lookup
STRATA = {
    "allmapp":      REGIONS / "GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz",
    "segdupmapp":   REGIONS / "GRCh38_notinlowmappabilitysegdups.bedsort.nochr.bed.gz",
    "nosegdupmapp": REGIONS / "GRCh38_notinlowmappabilitynosegdups.bedsort.nochr.bed.gz",
}
strata_starts, strata_ends = {}, {}
for s, fp in STRATA.items():
    ss, ee = [], []
    with gzip.open(fp, "rt") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) < 3 or cols[0] != CHR: continue
            ss.append(int(cols[1])); ee.append(int(cols[2]))
    order = sorted(range(len(ss)), key=lambda i: ss[i])
    strata_starts[s] = [ss[i] for i in order]
    strata_ends[s]   = [ee[i] for i in order]

def in_stratum(pos0, s):
    a, b = strata_starts[s], strata_ends[s]
    i = bisect.bisect_right(a, pos0) - 1
    return i >= 0 and pos0 < b[i]

# 3. AF bin assignment (matches __1_+__10_ stratification)
AF_BINS = [
    (0.001, 0.005, "0.001_0.005"),
    (0.005, 0.01,  "0.005_0.01"),
    (0.01,  0.05,  "0.01_0.05"),
    (0.05,  0.1,   "0.05_0.1"),
    (0.1,   0.5,   "0.1_0.5"),
    (0.5,   1.01,  "0.5_2"),
]
def af_bin(af):
    for lo, hi, lbl in AF_BINS:
        if lo <= af < hi: return lbl
    return None

# 4. Pre-load GA4K SNVs at AF >= 0.001
re_af = re.compile(r"AF=([0-9.eE\-]+)")
re_vt = re.compile(r"variant_type=([a-z]+)")
ga4k = {}
with gzip.open(GA4K_VCF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip().split("\t")
        if len(p) < 8: continue
        try: pos = int(p[1])
        except: continue
        ref, alt, info = p[3], p[4], p[7]
        if len(ref) != 1 or len(alt) != 1: continue
        vt = re_vt.search(info)
        if vt and vt.group(1) != "substitution": continue
        af_m = re_af.search(info)
        if not af_m: continue
        try: af = float(af_m.group(1))
        except: continue
        if af <= 0: continue
        ga4k[(pos, ref, alt)] = af
print(f"chr{CHR}: {len(ga4k)} GA4K SNVs preloaded", flush=True)

# 5. Build per-position gnomAD dicts (per-position counting, matches perl)
gnomad_pos_afbins  = collections.defaultdict(set)  # pos -> {af_bin} that has >=1 variant
gnomad_pos_alleles = collections.defaultdict(dict) # pos -> {(ref,alt): af}
n_gn = 0
with gzip.open(GNOMAD, "rt") as f:
    for line in f:
        cols = line.rstrip().split("\t")
        if len(cols) < 5: continue
        try: pos = int(cols[2])
        except: continue
        vid = cols[3]
        try: gn_af = float(cols[4])
        except: continue
        if gn_af <= 0: continue
        ab = af_bin(gn_af)
        if ab is None: continue
        try:
            ref_alt = vid.rsplit("_", 1)[-1]
            ref, alt = ref_alt.split("/")
        except: continue
        if len(ref) != 1 or len(alt) != 1: continue
        n_gn += 1
        gnomad_pos_afbins[pos].add(ab)
        gnomad_pos_alleles[pos][(ref, alt)] = gn_af

ga4k_pos_afbins = collections.defaultdict(set)
for (pos, ref, alt), ga_af in ga4k.items():
    ab = af_bin(ga_af)
    if ab is None: continue
    ga4k_pos_afbins[pos].add(ab)

print(f"chr{CHR}: scanned {n_gn} gnomAD SNVs over {len(gnomad_pos_afbins)} unique positions", flush=True)
print(f"chr{CHR}: {len(ga4k_pos_afbins)} GA4K positions with binnable AF", flush=True)

# 6. Tabulate per (stratum, k_class, af_bin) - per-position counting
ann_n    = collections.Counter()         # gnomAD positions
conf_n   = collections.Counter()         # positions with >=1 matched (ref,alt) in GA4K
lrs_n    = collections.Counter()         # GA4K positions
sum_logr = collections.defaultdict(float) # per-allele sum log10(gn/ga) for confirmed pairs
n_logr   = collections.Counter()         # count of allele pairs contributing

n_conf_total = 0
for pos, af_bins_at_pos in gnomad_pos_afbins.items():
    k_class = str(pos_to_k.get(pos, "nonGC"))
    pos0 = pos - 1
    strata_in = [s for s in STRATA if in_stratum(pos0, s)]
    if not strata_in: continue
    alleles_here = gnomad_pos_alleles[pos]
    for ab in af_bins_at_pos:
        confirmed_pos = False
        for (ref, alt), af_g in alleles_here.items():
            if af_bin(af_g) != ab: continue
            if (pos, ref, alt) in ga4k:
                confirmed_pos = True
                break
        for s in strata_in:
            key = (s, k_class, ab)
            ann_n[key] += 1
            if confirmed_pos:
                conf_n[key] += 1

for pos, alleles_here in gnomad_pos_alleles.items():
    k_class = str(pos_to_k.get(pos, "nonGC"))
    pos0 = pos - 1
    strata_in = [s for s in STRATA if in_stratum(pos0, s)]
    if not strata_in: continue
    for (ref, alt), af_g in alleles_here.items():
        ab = af_bin(af_g)
        if ab is None: continue
        ga_af = ga4k.get((pos, ref, alt))
        if ga_af is None: continue
        n_conf_total += 1
        lr = math.log10(af_g / ga_af)
        for s in strata_in:
            key = (s, k_class, ab)
            sum_logr[key] += lr
            n_logr[key] += 1

for pos, af_bins_at_pos in ga4k_pos_afbins.items():
    k_class = str(pos_to_k.get(pos, "nonGC"))
    pos0 = pos - 1
    strata_in = [s for s in STRATA if in_stratum(pos0, s)]
    if not strata_in: continue
    for ab in af_bins_at_pos:
        for s in strata_in:
            lrs_n[(s, k_class, ab)] += 1

print(f"chr{CHR}: confirmed gnomAD-GA4K allele matches = {n_conf_total}", flush=True)

# 7. Write per-cell counters
out_path = OUT_DIR / f"chr{CHR}.cell_counts.tsv"
with open(out_path, "w") as o:
    o.write("stratum\tk_class\taf_bin\tn_annotated\tn_confirmed\tsum_log_ratio\tn_lrs\tn_logr\n")
    all_keys = set(ann_n) | set(conf_n) | set(lrs_n) | set(n_logr)
    for key in sorted(all_keys):
        s, k, ab = key
        o.write(f"{s}\t{k}\t{ab}\t{ann_n[key]}\t{conf_n[key]}\t"
                f"{sum_logr[key]:.6f}\t{lrs_n[key]}\t{n_logr[key]}\n")
print(f"chr{CHR}: wrote {out_path}", flush=True)

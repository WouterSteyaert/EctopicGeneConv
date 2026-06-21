# 18 — Orthogonal long-read validation against GA4K PacBio HiFi

Orthogonal validation that the variant-density enrichment at gene-conversion-compatible positions is not driven by short-read paralog-mismapping artefacts in gnomAD. Reproduces the gnomAD-based enrichment formula and per-cell denominators from steps 06 + 07 (`__1_+__10_ContingencyStatistics`) using PacBio HiFi joint calls from the Genomic Answers for Kids cohort (n = 541 unrelated samples) as an independent variant source.

## Input data

| What | Where it must live | Source |
|---|---|---|
| GA4K joint HiFi VCFs (24 chrs, `pb_joint_merged.snv.chr*.vcf.gz` + tabix) | `$PROJECT_ROOT/geneconv_complete/artefact_validation_ga4k/vcf/` | <https://github.com/ChildrensMercyResearchInstitute/GA4K/tree/main/pacbio_snv_vcf> |
| gnomAD v4.1 genome BEDs (per chr, with AF column) | `$PROJECT_ROOT/geneconv_complete/gnomad/genome/{chr}.sorted.bed.gz` | Step 02 |
| `positions.gcdiffs` set-difference files per (k, chr) | `$PROJECT_ROOT/geneconv_complete/gcdiffs/{k}_{nextk}/{chr}.positions.gcdiffs.sorted.bed.gz` | Step 01 (`__7`–`__9`) |
| `conv/91/{chr}.repsum.geneconv.sorted.bed.gz` (longest-k file) | `$PROJECT_ROOT/geneconv_complete/conv/91/` | Step 01 |
| Mappability masks (GIAB v3.5: allmapp / segdupmapp / nosegdupmapp) | `$PROJECT_ROOT/geneconv_complete/regions/` | Step 00 setup |
| `_All_..._gnomad~genome_StatsSummary.R.out.txt` per stratum | `$PROJECT_ROOT/geneconv_complete/stats_{allmapp,segdupmapp,nosegdupmapp}/` | Steps 06 + 07 |

### Download GA4K VCFs

```bash
mkdir -p $PROJECT_ROOT/geneconv_complete/artefact_validation_ga4k/vcf
cd $PROJECT_ROOT/geneconv_complete/artefact_validation_ga4k/vcf
for chr in {1..22} X Y; do
    wget "https://raw.githubusercontent.com/ChildrensMercyResearchInstitute/GA4K/main/pacbio_snv_vcf/pb_joint_merged.snv.chr${chr}.vcf.gz"
    wget "https://raw.githubusercontent.com/ChildrensMercyResearchInstitute/GA4K/main/pacbio_snv_vcf/pb_joint_merged.snv.chr${chr}.vcf.gz.tbi"
done
```

Total ~540 MB. Files are sites-only VCFs with `INFO/AC`, `INFO/AN`, `INFO/AF`, `INFO/gnomAD_freq`, joint-called with DeepVariant + GLnexus on GRCh38.

## Pipeline

```
__1_PerChrVariantJoin.py             worker: per-chr (chr,pos,ref,alt) join of gnomAD x GA4K
__1_SubmitPerChrVariantJoin.pl       emits one SLURM job per chr (chr1-22, X, Y)
__2_AggregatePerCell.py              per-cell + per-AF contrast reference tables
__2_AggregateConfRates.py            pooled GC/non-GC confirmation rates (Wilson 95% CI)
__2_AggregateMhRd.py                 AF-adjusted Δ_mismatch + worst-case enrichment-inflation bound
__2_AggregateEnrichment.py           gnomAD vs LRS variant-density enrichment per (stratum, k, AF)
__3_MakeConfirmationFigure.R         headline figure (paired conf-rate lines)
__3_MakeEnrichmentFigure.R           enrichment lines + LRS/gnomAD ratio panel
__3_MakeForestFigure.R               forest plot Δ_mismatch + max-inflation panel
```

### Method

For every gnomAD genome SNV at AF ≥ 0.001, the worker assigns:
1. The set-difference longest-k class from `positions.gcdiffs` (k ∈ {17, 19, 21, 31, 41, 51, 61, 71, 81, 91} or `nonGC`).
2. The mappability stratum (same masks as the gnomAD contingency analysis).
3. The AF bin (6 bins matching `__1_+__10_ContingencyStatistics`).

Per (stratum, k_class, AF_bin) cell, counting is **per-position** (matches the perl `$Variants{$VarPos}` logic — a multi-allelic position with two alleles in the same AF bin counts once):
- `n_annotated` — gnomAD positions
- `n_confirmed` — positions where ≥ 1 (REF, ALT) matches the GA4K joint VCF
- `n_lrs` — GA4K positions (regardless of gnomAD overlap)
- `sum_log_ratio` — Σ log10(gnomAD_AF / GA4K_AF) summed over confirmed allele pairs (per-allele, for the AF-bias diagnostic)

The aggregators reconstruct the hierarchical paper partition (GC at k = K, non-GC at k = positions with longest_k < K + pure-nonGC) and apply the **identical formula and per-cell position denominators** as `__1_+__10_ContingencyStatistics`. As a result, `enr_gnomad` in `v2_enrichment_compare.tsv` reproduces the published gnomAD enrichment in `stats_<context>/_All_..._gnomad~genome_StatsSummary.R.out.txt` within rounding error (verified < 0.005 % over all 180 cells = 3 strata × 10 k × 6 AF).

`enr_lrs` is the same enrichment metric with GA4K SNVs as variant source — the orthogonal HiFi reproduction of the gnomAD enrichment.

### Run

```bash
export PROJECT_ROOT=/path/to/your/project

# 1. Emit per-chr SLURM jobs (24 chrs)
perl __1_SubmitPerChrVariantJoin.pl --ConfigFile=../00_Configuration/config.GRCh38.ini
for f in $PROJECT_ROOT/geneconv_complete/jobs/orthogonal_val_chr*.job; do sbatch $f; done

# 2. After all 24 jobs complete, aggregate (parallel, seconds each)
python3 __2_AggregatePerCell.py
python3 __2_AggregateConfRates.py
python3 __2_AggregateMhRd.py
python3 __2_AggregateEnrichment.py

# 3. Render figures (R 4.x + ggplot2, dplyr, tidyr, scales)
Rscript __3_MakeConfirmationFigure.R
Rscript __3_MakeEnrichmentFigure.R
Rscript __3_MakeForestFigure.R
```

### Output

All under `$PROJECT_ROOT/geneconv_complete/artefact_validation_ga4k/`.

```
results/
├── per_chr_v2/chr{1..22,X,Y}.cell_counts.tsv       per-chr worker output
├── v2_per_cell_aggregated.tsv                      per-cell raw counts
├── v2_contrast_per_af.tsv                          per-cell GC vs non-GC contrast
├── v2_conf_rates_reliableAF.tsv                    pooled conf rates (Wilson 95% CI) — for Confirmation figure
├── v2_mh_per_af.tsv                                per-AF Δ_mismatch + 95% CI
├── v2_mh_pooled.tsv                                AF-adjusted Δ_mismatch + max enrichment inflation per (stratum, k)
└── v2_enrichment_compare.tsv                       gnomAD vs LRS enrichment + ratio per (stratum, k, AF)

figures/
├── Confirmation_GC_vs_nonGC.{pdf,png}              headline figure
├── Forest_DeltaMismatch_MH.{pdf,png}               effect-size + 95% CI
├── MaxEnrichmentInflation.{pdf,png}                worst-case inflation bound (companion to forest)
├── Enrichment_gnomAD_vs_LRS.{pdf,png}              enrichment reproduction (paired lines per AF bin)
└── EnrichmentRatio_LRS_over_gnomAD.{pdf,png}       LRS/gnomAD ratio per cell
```

### Memory note

The per-chr worker pre-loads the per-chromosome positions.gcdiffs union into a dict; chr1 and chr2 (~180–195 M positions each) need ~24 GB RAM at peak. The submitter requests 32 GB to be safe; the rest fit easily in 12 GB.

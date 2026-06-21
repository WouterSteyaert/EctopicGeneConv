# Code for: Sequence architecture drives allele frequency through ectopic gene conversion in humans

## Overview

This repository contains all analysis scripts used in the study. Scripts are organized by analysis stage, matching the Methods section of the paper.

## Quick start

```bash
# 1. Point the pipeline at your local data directory:
export PROJECT_ROOT=/path/to/your/project        # see 00_Configuration/README.md for layout

# 2. Submit step 01 (k-mer extraction + GC-compatible position identification):
cd Methods_Scripts/01_IdentificationOfConversionCompatiblePositions
perl __1_SubmitFetchSequenceRepeatsFromFasta.pl --ConfigFile=../00_Configuration/config.GRCh38.ini
for f in $PROJECT_ROOT/geneconv_complete/jobs/chr*.rep.job; do sbatch $f; done
```

Every worker script reads paths and parameters from `00_Configuration/config.GRCh38.ini`
via the `ProjectConfig` helper. The placeholder `__PROJECT_ROOT__` in the config is
substituted at run-time with `$ENV{PROJECT_ROOT}` (or whatever is given through
`--ProjectRoot=...`). No file in the repository contains a hard-coded user path.

## Directory structure

| Directory | Description | Key scripts |
|-----------|-------------|-------------|
| 00_Configuration | Single source of truth for all paths and parameter lists used by step 01 (k-mer extraction). `config.GRCh38.ini` plus the `ProjectConfig.pm` helper that workers and submitters use to load it. | `config.GRCh38.ini`, `ProjectConfig.pm`, `README.md` |
| 01_IdentificationOfConversionCompatiblePositions | K-mer extraction, gene conversion-compatible position identification, set differences. Each worker `__N_*.pl` has a matching submitter `__N_Submit*.pl` that iterates over chromosomes × k-mer lengths × seeds and emits SLURM jobs. | Workers (`__1`–`__9`) + submitters (`__1_Submit`–`__9_Submit`) |
| 02_VariantData | gnomAD v4.1 genome variant extraction and formatting | Perl |
| 03_RandomVariantGeneration | Random control variant generation. One set of ~10 M uniformly sampled GRCh38 single-nucleotide substitutions is produced (directory named `random3` to match downstream analyses). Per-chrom counts scale by chromosome length (`<random>n_variants_total / genome_size_bp`); reference bases come from `<paths>reference_2bit` via `twoBitToFa`; the concat step restricts the combined set to the GIAB v3.5 all-mappable mask (`VariantsToQuery.bed`). Pass `--SetSeeds=` to the submitter for reproducible RNG. | Perl |
| 04_DeNovoMutationData | DNM extraction from Genomics England, quality filtering, Unfazed phasing. **Runs on the Genomics England Research Environment (LSF), not on the Ghent SLURM cluster.** Step order: `0_FindTriosToExclude` → `1_CreateListsOfDenovoVcfs` (emits LSF jobs for `2_FilterDNMs` per sublist) → `3_ProcessDNMs` → `4_UpdateDnmsWithGnomAd --Origin=<set>` → `5_SubmitPrepareInputUnfazed --Mode=fresh\|retry\|retry2` (one LSF job per trio, runs `5_PrepareInputUnfazed`) → `5_CheckAndSummarizeUnfazed`. Paths to the Genomics England LabKey TSVs and Research Environment are set via the `<paths>` keys `denovo_cohort_file`, `labkey_rd_analysis`, `labkey_rd_interpreted`, `denovo_work_dir`, and `reference_fasta_full` in `config.GRCh38.ini`. | Perl |
| 05_GWASVariantData | NHGRI-EBI GWAS Catalog ingest. The submitter splits the full catalog TSV into N shards and emits one SLURM job per shard (`__5_FetchGwasVariants`), each resolving the catalog's rsIDs to GRCh38 coordinates via the Ensembl variation API (release set by `<gwas>ensembl_release`) and classifying each catalog risk allele as Risk or Protection. `__5_WriteCleanGwasVariants` then merges all shards, joins each variant to its gnomAD AF (per-chromosome BED from step 02), and writes per-chromosome BEDs stratified by `<class>_<type>` (Disease_Risk, Disease_Protection, Trait_All) plus a combined `All/` set; the phenotype→class map is read from `classified_phenotypes`. Requires the public Ensembl Perl API. | Perl |
| 06_EnrichmentAnalysis | 2×2 contingency-statistics enrichment per (k × AF bin × mappability context) combination. `__1_ComputeContingencyStatistics` is the per-window worker (selects context via `--MappabilityContext=allmapp\|segdupmapp\|nosegdupmapp` and variant set via `--QueryVarList` + `--AfColumn`); `__1_SubmitComputeContingencyStatistics` iterates chr × k × AF-bin × context and emits one SLURM job per cell (defaults to the gnomAD genome case; override `--QueryVarList`/`--AfColumn`/`--VarFreqMode=unstratified` for DNM, GWAS, random sets). `__10_SummarizeContingencyStatistics` aggregates the per-window worker outputs across all chromosomes and writes a three-block `*.pos.txt` summary (PercDiff_Pos, ExcessPercPos, ExcessPercVar) plus an R driver that performs the chi-squared, binomial, and concordance tests. | Perl + R |
| 07_MultipleTestingCorrection | Bonferroni correction on the per-(k, AF-bin) p-values produced by step 06. `CorrectForMht.pl` reads every `*.R.out.txt` in each `<StatsRoot>/stats_<context>/` directory and emits an R driver that sets `m = nrow(Table)` (90 for stratified gnomAD / DNM / GWAS; 10 for unstratified random) and applies `p.adjust(method="bonferroni")` to both the positional and concordance chi-squared p-values, plus a log10(m) shift for the log-scale exponents. Each variant set is corrected independently. Pass `--StatsRoot=<dir>` to point at a non-default stats tree (e.g. the DNM export root). | Perl + R |
| 08_SlidingWindowAnalysis | Sliding-window enrichment (1 Mb / 50 kb step) on the step-06 contingency outputs. `1_correct_enrichment.pl` is the paper version: per-window odds ratio with pseudocount 0.5, control group corrected by subtracting GC-compatible counts from all longer k. `2_identify_peaks.pl` finds non-overlapping 1 Mb outlier windows (mean + 2 SD on log2 OR) and merges adjacent outliers into peaks (≥ 10 gc_var, annotated with SD overlap + pericentromeric/subtelomeric/interstitial classification). `3_excl_peaks_table.sh` rebuilds the genome-wide enrichment table after excluding peak windows (Supplementary Table 12). `__1_ExtractSlidingWindowEnrichment.pl` is a diagnostic uncorrected variant retained for comparison. | Perl + bash |
| 09_RecombinationHotspotAnalysis | Enrichment of variants at conversion-compatible positions inside vs outside meiotic recombination hotspots. Flow: `__1_PrepareHotspotBED.sh` converts the deCODE sex-averaged bigWig (`recombAvg.bw`) to bedGraph and thresholds at ≥10 cM/Mb to produce the hotspot BED; `__7_SubmitContingencyStatisticsHotspot.pl` emits per (chrom × k × AF-bin × chunk) SLURM jobs that run `__7_ComputeContingencyStatisticsHotspot.pl` (step-06 worker logic restricted to hotspot-overlapping positions); `__10_MergeChunkedRecombHotspots.pl` merges the chunked outputs; `SupTable_RecombHotspots.R` computes log₂ OR inside vs outside, the Woolf-style z-test for heterogeneity, and Bonferroni correction (m = 90). | Perl + R |
| 10_DistanceAnalysis | Donor-acceptor distance analysis. Flow (run per k via `--REPLEN=`): `1_filter_1to1_coords.sh` (keep 1:1 pairs from the conversion coords files) → `2_extract_1to1_distances.sh` (dedup pairs by canonical ordering, classify intra/inter, distance bins) → `3_compute_pair_gcdiffs_v2.sh` (set-difference: keep pairs exclusive to this k) → `4_annotate_gnomad_diff.sh` (annotate both positions with allele-matched gnomAD AF) → `5_summary_tables.sh` (per-(k × distance) summaries) → `6_plots.R` + `Figure4_DistanceFunctional.R` (figure; includes the leave-one-k-out CV linear regression of mean AF on k and log₁₀(distance midpoint)). All paths derive from `$PROJECT_ROOT`. The GC-biased-conversion analysis built on these per-pair outputs lives in step 11. | bash + R |
| 11_GCBiasAnalysis | GC-biased gene conversion analysis on the step-10 per-pair gnomAD-annotated outputs. `7_gbgc_summary.sh` classifies each variant as W→S (A/T → G/C), S→W (G/C → A/T) or other and writes three tables: (i) mean/median AF per (k × direction), (ii) mean AF per (k × distance bin × direction), and (iii) fraction of W→S among all W↔S variants at cumulative AF thresholds (0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5). `7_gbgc_plots.R` renders the figures. | bash + R |
| 12_LDAnalysis | Linkage disequilibrium comparison restricted to 1:1 paralog pairs. Flow: `__8_Extract1to1Positions` (extract 1:1 GC positions per k from step-10 `distances.diff` output) → `__9_ClassifyVariantsForLd_1to1` (classify gnomAD variants into GC/non-GC per mappability stratum × k × MAF bin; non-GC excludes 1:1 GC positions from ALL k) → `__10_MergeChromosomes_1to1` → `__11_SubmitLdCalculation_1to1` (per (mapp × pop × k × MAF bin): PLINK `--snps-only just-acgt`; per-focal mean r² across partners within 10/25/50/100 kb; GC focals exclude partners within ±(k−1)/2 bp because the flanking k-mer is identical between donor and acceptor by construction; max 5,000 variants per chr per class). `__12_AnalyzeLd_1to1.R` performs the Wilcoxon rank-sum test between GC and frequency-matched non-GC variants. `__13_LD_vs_RecombRate_1to1.R` stratifies by deCODE recombination rate (six bins). EUR (503 samples) and ALL (3,202 samples) populations are analysed separately. | Perl + R |
| 13_TractLengthValidation | BLAST-based homologous tract detection and validation. `__2_SubmitComputeTractLengthWgs` samples N gnomAD variants per AF bin (`<tract>n_variants_per_bin`, default 5000 × 9 bins = 45,000) from `gnomad_genome_dir/VariantsToQuery.bed` and batches them into SLURM jobs that call `__2_ComputeTractLengthWgs`; `__1_SubmitComputeTractLengthRandom` does the same for random control positions (`<tract>n_random_variants`, default 5000). Each worker iterates over `<tract>window_sizes` (default 20-90 bp half-windows), introduces the alternative allele in the centre, BLASTs against `<paths>reference_fasta` with `blastn -dust no -word_size <tract>blast_word_size -perc_identity <tract>blast_perc_identity`, parses each alignment to identify flanking single-unique-nucleotides (SUNs), and emits per-alignment tract length + tract analysis length. `__3_SubtractRandomFromWgs` builds the AF-binned tract-length distribution after subtracting the random baseline; `__4_CompareTractLengthWithGcDiff` joins BLAST tract analysis lengths to the k-mer-based gcdiff classification per chromosome. | Perl |
| 14_PerGenerationRate | Per-generation ectopic conversion rate (r<sub>gc</sub>) estimation from DNM contingency tables (step 06 + step 07). For each (k × mappability stratum), computes the excess of DNMs at conversion-compatible positions over the discordant background, then r<sub>gc</sub> = Excess / (N × ConcorTotal), with N = 11,963 trios. Two CIs are emitted: (1) delta-method on log(fold) for fold excess and r<sub>gc</sub>/μ<sub>eff</sub>; (2) Poisson Var(Excess) = ConcorVar + (ConcorTotal/DiscorTotal)² × DiscorVar for r<sub>gc</sub> itself. Negative lower bounds are clipped to 0 (see Methods for the k=17 SD edge case). Reads from `<denovo_work_dir>/export/stats_<context>/`; outputs `per_generation_gc_rate.tsv` and a two-panel PDF. | R |
| 15_ComputationalImplementation | Software versions and compute environment description | — |
| 16_Figures | Main figure generation scripts: `Figure2_EnrichmentDNM.R` (Fig 2 — DNM and variant enrichment heatmaps), `Figure3.R` (Fig 3 — sliding-window genome-wide signal), `Figure5_LD_GWAS.R` (Fig 5 — LD comparison + GWAS depletion). All read from `$PROJECT_ROOT/geneconv_complete/...` and write PDFs to `$PROJECT_ROOT/geneconv_complete/figures/`. Figure 4 (donor-acceptor distance) lives next to its analysis pipeline in step 10 (`Figure4_DistanceFunctional.R`); the per-generation r<sub>gc</sub> figure lives in step 14. | R |
| 17_SupplementaryTables | Supplementary table generation. `SupTable_SlidingWindowPeaks.R` produces ST11 from the corrected peak set (step 08, MIN_GC_VAR = 10). Other supplementary tables are produced next to their analyses: ST7-9 from step 06, ST10 from step 14, ST12 from step 08, ST13 from step 09, ST14 from step 10, ST16/ST18 from step 12. The reversal-regime supplementary tables (former ST19) were removed after the LD methodological correction. | R |
| 18_OrthogonalLongReadValidation | Orthogonal validation that the variant-density enrichment is not driven by short-read paralog mismapping in gnomAD. Reproduces the gnomAD enrichment formula and per-cell denominators from steps 06 + 07 using PacBio HiFi joint calls from Genomic Answers for Kids (n = 541 unrelated samples; <https://github.com/ChildrensMercyResearchInstitute/GA4K>) as an independent variant source. `__1_PerChrVariantJoin.py` (per-chr worker) joins gnomAD x GA4K at the (chr, pos, REF, ALT) level and tabulates per-position counts per (mappability stratum, longest-k class, AF bin). `__2_AggregateEnrichment.py` reconstructs the hierarchical paper partition (matches `__1_+__10_ContingencyStatistics` exactly: enr_gnomad reproduces published values within < 0.005 % over all 180 cells) and computes the orthogonal HiFi-based enrichment `enr_lrs`. `__2_AggregateMhRd.py` computes the AF-adjusted Delta_mismatch and worst-case enrichment-inflation bound. `__3_Make*Figure.R` render the headline confirmation figure, enrichment-reproduction lines, and forest plot. | Python + Perl + R |

## Dependencies

- **Perl** v5.34.0 with modules: Config::General, IO::Zlib, Getopt::Long
- **R** v4.3.2 with packages: ggplot2, dplyr, tidyr, cowplot, pwr
- **BEDTools** v2.30.0
- **htslib/tabix** v1.15
- **BCFtools** v1.18
- **PLINK** v1.9
- **BLAST+** (blastn)
- **twoBitToFa** (UCSC utilities)
- **Unfazed** (DNM parent-of-origin phasing)
- **Ensembl Perl API** (release matching `<gwas>ensembl_release`; required by step 05 for rsID-to-coordinate resolution; connects to the public Ensembl MySQL server)

## Input data

All public datasets used in the paper, with the exact download URLs verified to
return HTTP 200 at the time of writing. Place each dataset under the path expected
by `config.GRCh38.ini` (see the `<paths>` section).

### GRCh38 reference genome
- `hg38.2bit`: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit>
- `hg38.fa.gz`: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>
  (gunzip; downstream scripts also expect per-chromosome FASTA `chr{N}.fa.gz` —
  split with `samtools faidx hg38.fa chr{N} | bgzip > chr{N}.fa.gz`)

### gnomAD v4.1 genomes (per-chromosome sites VCF)
- Pattern: `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{N}.vcf.bgz`
- Example: `gnomad.genomes.v4.1.sites.chr1.vcf.bgz`, ..., `chr22.vcf.bgz`, `chrX.vcf.bgz`, `chrY.vcf.bgz`

### 1000 Genomes high-coverage 30x NYGC, phased
- Directory: <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/>
- Per-chromosome pattern: `1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`

### GIAB genome stratifications v3.5 (mappability + segdup BEDs)
- Archive: <https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.5/genome-stratifications-GRCh38@all.tar.gz>

### NHGRI-EBI GWAS Catalog, release 2024-11-20 (Ensembl 113)
- File: <http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2024/11/20/gwas-catalog-associations.tsv>
- The paper text refers to this file as `gwas_catalog_v1.0.2-associations_e113_r2024-11-20.tsv` (renamed locally to record the version + Ensembl release in the filename); the bytes are identical.

### deCODE sex-averaged recombination map
- bigWig: <https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw>
- Step 09 (`__1_PrepareHotspotBED.sh`) converts to bedGraph with `bigWigToBedGraph`.

### Genomic Answers for Kids (GA4K) PacBio HiFi joint VCFs (step 18)
- Per-chromosome pattern: `https://raw.githubusercontent.com/ChildrensMercyResearchInstitute/GA4K/main/pacbio_snv_vcf/pb_joint_merged.snv.chr{N}.vcf.gz`
- 24 chromosomes (chr1-22, X, Y), ~540 MB total; sites-only, INFO/AF + INFO/gnomAD_freq, joint-called with DeepVariant + GLnexus (n = 541 unrelated samples; GRCh38).
- Place under `$PROJECT_ROOT/geneconv_complete/artefact_validation_ga4k/vcf/`. See `18_OrthogonalLongReadValidation/README.md` for the download loop.

### Genomics England 100,000 Genomes Project DNMs (restricted access)
- Available only inside the Genomics England Research Environment via LabKey;
  paths set under `<paths>denovo_*` and `<paths>labkey_*` in `config.GRCh38.ini`.
- The pipeline step that uses these files (`04_DeNovoMutationData/`) runs on the
  GEL LSF cluster, not on the public-data side.

## Execution order

Scripts within each directory are numbered and should be executed in order. Dependencies between directories follow the numbering: later directories may depend on output from earlier ones. The enrichment analysis (06) requires output from both the position identification (01) and variant data (02-05). The sliding window (08), distance (10), LD (12), and other downstream analyses require output from the enrichment framework.

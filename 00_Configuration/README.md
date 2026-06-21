# Configuration

## Quick start

```bash
# 1. Set PROJECT_ROOT to the directory containing your reference + working data
export PROJECT_ROOT=/path/to/your/project

# 2. Submit scripts read the config:
perl 01_IdentificationOfConversionCompatiblePositions/__1_SubmitFetchSequenceRepeatsFromFasta.pl \
     --ConfigFile=00_Configuration/config.GRCh38.ini
```

## What the config holds

- **paths**: project root, reference FASTA pattern (per chr), reference dict,
  and data subdirectories (`repeats/`, `repsum/`, `conv/`, `gcdiffs/`, `jobs/`,
  `regions/`).
- **genome**: assembly label, chromosome list, chunk size.
- **repeat_lengths**: which k-mer lengths the pipeline scans, and which short-long
  pairs are used by the set-difference step (`__7`).
- **slurm**: default per-job SLURM resources.
- **modules**: software modules to `module load` inside each worker job.

The token `__PROJECT_ROOT__` in `config.GRCh38.ini` is replaced at submit-time
by `$ENV{PROJECT_ROOT}`. If you copy the config to a new project, leave the
token in place — only `PROJECT_ROOT` needs to change.

## External input data

Place the following under `$PROJECT_ROOT/`:

| Path | Source | Notes |
|------|--------|-------|
| `hg38/chr1.fa.gz` ... `chr22.fa.gz`, `hg38/hg38.dict` | UCSC Genome Browser (GRCh38) | per-chromosome gzipped FASTA + Picard-style sequence dictionary |
| `regions/GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz` (+ `.tbi`) | GIAB v3.5 genome stratification | bgzip + tabix-indexed; chromosome labels without `chr` prefix |
| `regions/GRCh38_notinlowmappabilitysegdups.bedsort.nochr.bed.gz` | GIAB v3.5 | as above (segdup-only mappable) |
| `regions/GRCh38_notinlowmappabilitynosegdups.bedsort.nochr.bed.gz` | GIAB v3.5 | as above (non-segdup mappable) |
| `extra/Variation/gnomAD/genomes/gnomad.genomes.v4.1.sites.chrN.vcf.bgz` (+ `.tbi`) | gnomAD v4.1 release | per-chromosome PASS-filtered VCFs; input for step 02 |

Preprocessing the GIAB BEDs:

```bash
# Download (example, see GIAB v3.5 release notes for exact URL)
wget https://.../GRCh38_notinlowmappabilityall.bed.gz
zcat GRCh38_notinlowmappabilityall.bed.gz | sed 's/^chr//' \
    | sort -k1,1 -k2,2n | bgzip > GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz
tabix -p bed GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz
```

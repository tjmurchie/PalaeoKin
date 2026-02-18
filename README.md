# PalaeoKin

**PalaeoKin** is a lightweight, ancient-DNA-friendly pipeline for estimating **pairwise relatedness** from **mapped BAM files** at low and uneven coverage, with built-in **QC plots** that make it clear when/why estimates are unstable.

It wraps:
- **ANGSD** (genotype likelihoods + site filtering)
- **NgsRelate** (relatedness: KING-robust kinship, R0, R1)
- **samtools** (depth / diagnostic summaries)

Outputs include a human-readable **summary report** plus Illustrator-friendly **vector plots (PDF/SVG)**, with optional **downsampling** and **mixture/allele-balance** diagnostics.

> Version: **1.0.0**


## When to use PalaeoKin

PalaeoKin is useful when you:
- have **ancient/historical DNA** with short fragments and damage patterns
- have **low or uneven endogenous coverage**
- want a **BAM-in → report-out** workflow
- need QC that supports a defensible “**inconclusive**” call when data are too sparse

It’s particularly handy outside human SNP-panel workflows (e.g., wildlife archaeology / palaeogenomics).


## Quickstart

1) Make a BAM list (one BAM path per line):

```bash
printf "%s\n" /path/to/sample1.bam /path/to/sample2.bam > bam_list.txt
```

2) Run:

```bash
./PalaeoKin.sh \
  --bams bam_list.txt \
  --ref /path/to/reference.fa \
  --outdir ./palaeokin_out \
  --threads 8 \
  --downsample \
  --mixture-check
```

3) Key outputs:
- `palaeokin_out/reports/summary_report.txt`
- `palaeokin_out/plots/king_r0_plot.(png|pdf|svg)`
- `palaeokin_out/plots/downsampling_plot.(png|pdf|svg)`
- `palaeokin_out/plots/qc_panels.(png|pdf|svg)` (Panels A–D)


## Installation

### Dependencies
Tools expected on `PATH`:
- `samtools` (required)
- `angsd` (required)
- `NgsRelate` (required)

Python:
- Python 3 + `numpy`, `pandas`, `matplotlib`

> Note: ANGSD/NgsRelate installation varies by HPC. Many clusters provide them via modules.

### Conda

PalaeoKin will also attempt to create/verify the `palaeokin` environment automatically on first run via `scripts/00_setup_environment.sh`.


```bash
conda env create -f environment.yml
conda activate palaeokin
```

Then confirm tools:
```bash
samtools --version
angsd -h | head
NgsRelate | head
python -c "import numpy,pandas,matplotlib"
```


## Output structure

```
outdir/
  reports/
    summary_report.txt
    run_summary.json
  plots/
    king_r0_plot.png|pdf|svg
    downsampling_plot.png|pdf|svg
    qc_panels.png|pdf|svg     # A–D
  ngsrelate/
    relatedness_results.tsv
  validation/
    downsampling/
      downsampling_results.tsv
    mixture/
      <sample>/
        maf_values.tsv.gz
        maf_histogram.png|pdf|svg
        mixture_report.txt
    depth_at_informative_sites.tsv
  angsd/
    <param_set>/
      sites.txt
      joint_gl.glf.gz
      ...
```


## Plot formats (PNG + PDF + SVG)

By default, plots are written as **png,pdf,svg**.

To change formats:
```bash
export PLOT_FORMATS="png,pdf"
```

SVGs are post-processed to use **Arial** (to reduce Illustrator font warnings) while keeping text editable.


## Interpreting results

### Relatedness (KING / R0 / R1)
PalaeoKin reports:
- **KING-robust kinship** (interpretation, approximate):
  - ~0.00 = unrelated
  - ~0.0625 = 3rd degree
  - ~0.125 = 2nd degree
  - ~0.25 = 1st degree
  - ~0.50 = identical/duplicate

Thresholds used in plots:
- Identical ≥ 0.354
- 1st degree ≥ 0.177
- 2nd degree ≥ 0.0884
- 3rd degree ≥ 0.0442

It also reports:
- **R0**: fraction of homozygous-opposite sites (closer relatives generally have lower R0)
- **R1**: allele-sharing statistic used with R0 + KING

**Sites used matters most at low coverage. Rule of thumb:**
- `< ~100 sites`: exploratory only
- `< ~500 sites`: often unstable
- `> ~2000 sites`: typically much more reliable (assuming decent mapping/reference)

### QC panels (A–D)
- **A:** number of informative sites under different ANGSD parameter sets
- **B:** depth distribution at informative sites (shows which sample is limiting)
- **C:** minor allele fraction histogram (allele-balance heuristic)
- **D:** minor allele fraction vs depth (distinguishes low-depth discreteness vs mixture-like behaviour)

### Mixture / allele-balance (counts-based)
For each informative site:
\[
  \mathrm{MAF}=\frac{\mathrm{minor}}{\mathrm{major}+\mathrm{minor}}
\]
- MAF ranges **0–0.5** by definition (minor ≤ major).
- For a **single diploid individual**, heterozygous-like sites tend to cluster near **~0.5** at sufficient depth.
- At low depth (3–6×), discrete outcomes (1/3≈0.33, 2/5=0.40, etc.) are expected from sampling noise.

Mixture is more plausible if intermediate MAF values persist at **higher depth** (use Panel D).


## Tips / pitfalls
- **Reference choice matters:** cross-species mapping can reduce usable sites and introduce bias.
- **Coverage imbalance matters:** if one sample is extremely low coverage, estimates will be driven by that sample and may be unstable.
- **Transversions-only** is recommended for ancient DNA to reduce damage-driven artefacts.
- QC plots are designed to support defensible “inconclusive” calls.


## Citation
Please cite:
- ANGSD
- NgsRelate
- and this repository (see `CITATION.cff`)


## License
MIT (see `LICENSE`).

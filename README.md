# sSNV Calling

These Python scripts filter candidate somatic single-nucleotide variants (SNVs)
called against a phased, donor-specific assembly (DSA; contigs named
`haplotype1-*` / `haplotype2-*`). They combine evidence from PacBio HiFi and
Illumina reads across many tissues of one donor. The pipeline:

* removes calls in assembly regions flagged as unreliable (NucFlag, distortion maps)
* removes clustered calls and collapses duplicate calls from homologous or paralogous sequence
* removes variants common in gnomAD v4.1
* re-checks read support for each variant on each haplotype, per platform and per tissue
* merges the evidence across platforms and tissues, phases each variant, and applies
  binomial (VAF) and Fisher strand-bias tests with Benjamini–Hochberg correction
* removes variants in GA-rich or homopolymer sequence contexts

The full algorithm is described in the Methods section of the manuscript.

---

## Contents

```
README.md               this file
LICENSE                 MIT license
requirements.txt        Python dependencies
*.py                    pipeline source code (Python 3)
helpers/                shell wrappers used to run each step in production (HPC paths)
demo/data/              small synthetic demo dataset
demo/run_demo.sh        runs the pipeline on the demo dataset
demo/expected_output/   expected demo outputs
```

In the submission zip the Python scripts are in `scripts/` and the wrappers in
`scripts/helpers/`. `demo/run_demo.sh` works with either layout.

---

## 1. System requirements

### Software dependencies

| Dependency | Version | Needed by |
|---|---|---|
| Python | ≥ 3.9 (tested 3.12.13 and 3.14.3) | all scripts |
| numpy | ≥ 1.24 (tested 2.5.4) | scipy / statsmodels |
| scipy | ≥ 1.10 (tested 1.18.1) | `combine_tissues.py`, `single_tissue_combine_platform.py` |
| statsmodels | ≥ 0.14 (tested 0.15.0) | `combine_tissues.py`, `single_tissue_combine_platform.py` |
| tqdm | ≥ 4.60 (tested 4.70.1) | most scripts (progress bars) |
| pysam | ≥ 0.21 (tested 0.24.1) | `check_remapped_reads.py`, `collapse_duplicates.py`, `ga_richness_filter.py` |
| pybedtools + **bedtools** binary | pybedtools ≥ 0.9 (tested 0.12.1); bedtools ≥ 2.30 | `remove_flagged_regions.py`, `remove_distorted_regions.py`, `combine_callsets.py` |
| pandas | ≥ 1.5 (tested 3.0.6) | `combine_callsets.py` (via pybedtools) |
| hail, gnomad_toolbox, Java 11+ | latest | `annotate_gnomad.py` only (downloads gnomAD v4.1 genomes; needs network access) |
| bash, awk | any | `helpers/*.sh`, `demo/run_demo.sh` |

### Operating systems

* Linux or macOS. Windows is not supported (pysam and bedtools do not support it natively).
* **Tested on:**
  * Rocky Linux 8.10 (Green Obsidian), Python 3.12.13: production runs on an HPC cluster.
  * macOS 26.3 (Apple M1 Pro, 16 GB RAM), Python 3.14.3: the full demo.

### Hardware

No non-standard hardware is required. The demo runs on a normal desktop or
laptop. On full-size data, `check_remapped_reads.py` and `collapse_duplicates.py`
need enough RAM to hold the per-read metrics for all candidate sites (we used an
HPC node). `annotate_gnomad.py` starts a local Hail/Spark session with 16 cores.

---

## 2. Installation guide

```bash
# 1. Get the code
git clone https://github.com/mdnoyes/smart_p25.git && cd smart_p25

# 2. Create an isolated environment and install the Python dependencies
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 3. (Optional) the region filters need the bedtools binary
#    macOS:  brew install bedtools
#    Linux:  sudo apt-get install bedtools   or   conda install -c bioconda bedtools

# 4. (Optional) gnomAD annotation only
pip install hail git+https://github.com/broadinstitute/gnomad-toolbox
```

You can use conda instead of pip:
`conda create -n snvfilter -c conda-forge -c bioconda python=3.12 numpy scipy statsmodels tqdm pysam pybedtools pandas bedtools`.

The scripts are run directly with `python <name>.py`, so nothing needs to be compiled.

**Typical install time:** about 1–2 minutes on a normal desktop with a broadband
connection (binary wheels for every dependency). Installing Hail adds a few more minutes.

---

## 3. Demo

The demo runs the pipeline steps that do not need multi-gigabyte inputs (CRAMs,
full assemblies, gnomAD) on a small synthetic dataset in `demo/data/`:

| Stage | Scripts | Demo input |
|---|---|---|
| A. Clustered-variant filter | `inter_snp_distance.py` → `get_clusters.py` → `remove_clusters.py` | 12 candidate SNVs |
| B. gnomAD frequency filter | `filter_gnomad.py` | 6 gnomAD-annotated SNVs |
| C1. Per-tissue platform merge | `single_tissue_combine_platform.py` | HiFi + Illumina read counts for 8 SNVs × 3 tissues |
| C2. Cross-tissue merge and statistical tests | `combine_tissues.py` | output of C1 |
| C3. Test filter | `awk` (same logic as `helpers/remove_tests.sh`) | output of C2 |
| C4. GA-richness / homopolymer filter | `ga_richness_filter.py` | output of C3 + 2-contig, 4 kb demo DSA |
| D. Region filters (optional, needs bedtools) | `remove_flagged_regions.py`, `remove_distorted_regions.py` | 12 candidate SNVs + BED files |

Each synthetic variant in stage C is built to exercise one filter
(`v1_somatic_h1`, `v2_germline_like`, `v3_strand_bias`, `v4_low_confidence`,
`v5_recurrent_err`, `v6_ga_rich`, `v7_homopolymer`, `v8_somatic_h2`).

### Running the demo

```bash
bash demo/run_demo.sh
# or with a specific interpreter:
PYTHON=.venv/bin/python bash demo/run_demo.sh
```

Outputs go to `demo/output/`. The script compares each one with
`demo/expected_output/` and exits non-zero if any differ. Stage D is skipped
(reported as `SKIPPED`) if bedtools is not installed.

### Expected output

The script ends with:

```
=== Comparing outputs with demo/expected_output/
OK       A.clusters.tsv
OK       A.intersnp_annotated.tsv
OK       A.no_clustered_variants.tsv
OK       B.no_gnomad.tsv
OK       C1.blood.tsv
...
OK       C4.no_homopolymer.tsv
OK       D.no_distortions.tsv          (or SKIPPED without bedtools)
OK       D.no_flagged_regions.tsv      (or SKIPPED without bedtools)

Demo finished successfully. Final call set: .../demo/output/C4.no_homopolymer.tsv
```

What each stage should do:

* **A.** The four variants within 500 bp at chr1:1000–1450 are removed. The pair
  at chr1:8000/8300 and the three variants at chr2:3000–3300 are kept and annotated
  `in_cluster=True`. 8 of 12 variants remain.
* **B.** Variants with gnomAD AF ≥ 1e-4 (`g3`, `g4`) are removed. 4 of 6 remain.
* **C2.** `v4` (low confidence only) and `v5` (recurrent error on both haplotypes)
  are not reported. The 6 remaining variants are phased and tested.
* **C3.** `v2` (VAF ≈ 0.5, fails both binomial tests) and `v3` (Illumina strand
  bias) are removed.
* **C4.** `v6` (GA-rich context, HiFi-only support in one tissue) and `v7`
  (adjacent to an `AAAAA` homopolymer) are removed. **Final call set: `v1_somatic_h1`
  (h1; blood, liver) and `v8_somatic_h2` (h2; liver, lung).**
* **D.** NucFlag `collapse`/`misjoin` regions remove chr1:1200 and chr2:9000.
  chr1:5000 is kept with `NucFlag=het`. The distortion maps remove chr1:8000 and chr2:3100.

### Expected run time

About 5 seconds on a normal desktop (about 25 seconds on the first run while
Python compiles and caches the scipy/statsmodels imports).

---

## 4. Instructions for use

### Pipeline order

The production pipeline runs once per donor (`$sample`, e.g. `SMHT004`). The
wrappers in `helpers/` show the exact invocations we used; edit the
hard-coded HPC paths to match your system.

| # | Step | Script(s) | Main inputs | Output |
|---|---|---|---|---|
| 1 | Merge fibroblast and merged-tissue call sets | `combine_callsets.py -f FIBRO.vcf.gz -t TISSUE.vcf.gz -o OUT` | bgzipped VCFs (FILTER `PASS`/`GERMLINE`; FORMAT `GQ:DP:AD:VAF`) | candidate TSV |
| 2 | Remove NucFlag-flagged regions | `remove_flagged_regions.py -i IN -b NUCFLAG.bed -o OUT` | TSV with `chrom`,`pos`; NucFlag BED (4th column = flag) | `*.in_good_regions.tsv` |
| 3 | Remove clustered SNVs | `helpers/remove_clusters.sh` (`inter_snp_distance.py` → `get_clusters.py` → `remove_clusters.py`) | sorted TSV with `chrom`,`pos`,`id` (`id` = `chrom_pos_ref_alt`) | `*.no_clustered_variants.tsv` |
| 4 | Collapse duplicate calls (homolog liftover + alternate alignments) | `collapse_duplicates.py -v IN -b HOMOLOG.bed -f DSA.fasta -c FIBRO.cram -o OUT` | homolog BED (4th column = variant id), DSA FASTA (indexed), CRAM with `XA` tags | `*.collapsed_duplicates.tsv` |
| 5 | Remove distorted regions | `remove_distorted_regions.py -i IN -b1 HAP1.bed -b2 HAP2.bed -o OUT` | distortion-map BEDs | `*.no_distortions.tsv` |
| 6 | gnomAD annotation and filter | `annotate_gnomad.py -v IN -b HG38_LIFTOVER.bed -o OUT`, then `filter_gnomad.py -v IN -o OUT` | BED of GRCh38 positions (4th column = variant id) | `*.no_gnomad.tsv` |
| 7 | Re-check read support per tissue and platform | `check_remapped_reads.py -v IN -r READS.tsv.gz -p {hifi,illumina} -c TISSUE.cram -f DSA.fasta -o OUT` | candidate TSV (`variant_id`, `chrom`, `pos`, `ref`, `alt`, `homolog_*`, `collapsed_calls`); gzipped per-read metrics; CRAM with `AS`, `oh` and `XA` tags | per-tissue, per-platform read counts |
| 8 | Merge platforms within each tissue | `single_tissue_combine_platform.py -v VARIANTS -hifi HIFI_COUNTS -ilm ILM_COUNTS -o OUT` | output of step 7 (either platform may be omitted) | `tissue_data/$sample.$tissue.tsv` |
| 9 | Merge tissues; phasing; binomial and strand-bias tests | `combine_tissues.py -v VARIANTS -t NAME FILE [-t NAME FILE ...] -o OUT` | output of step 8 for every tissue | `*.combined_tissues.tsv` |
| 10 | Remove variants failing the tests | `helpers/remove_tests.sh $sample` | output of step 9 | `*.passed_strand_bias.tsv` |
| 11 | GA-richness / homopolymer filter | `ga_richness_filter.py -v IN -f DSA.fasta -o OUT` | output of step 10 | `*.no_homopolymer.tsv` (final) |

Run any script with `-h` to list its arguments.

### Input format notes

* All intermediate files are tab-separated with a single header row.
* Contig names must contain `haplotype1` or `haplotype2`; the haplotype is parsed
  from the first character after `haplotype`. Contigs containing `#mt#` are skipped
  in step 7.
* The per-read metrics file for step 7 (gzipped TSV) needs these columns:
  `read_name`, `chromosome`, `snv_position`, `read_length`, `base_quality`,
  `gap_compressed_identity`, `aligned_fraction`, `dist_to_read_end`, `base_at_snv`, `strand`.
* CRAMs for step 7 need the alignment score (`AS`), original haplotype (`oh`) and
  optional alternate alignments (`XA`) tags.
* Platform-specific read-quality and read-support thresholds are set in
  `ReadFilter.for_type` and `PresenceFilter.for_type` in `check_remapped_reads.py`.
* Tissue names passed to `combine_tissues.py -t` must not contain underscores
  (use hyphens, e.g. `brain-frontal`).

### Reproduction instructions

The manuscript call sets cover five donors: `SMHT001`, `SMHT004`, `SMHT005`,
`SMHT012` and `SMHT024`. For each donor, run the steps in the table above in
order. Steps that have a wrapper in `helpers/` can be run as
`bash helpers/<wrapper>.sh $sample` after editing the paths at the top of the wrapper:

```bash
for sample in SMHT001 SMHT004 SMHT005 SMHT012 SMHT024; do
    # steps 1-2: combine_callsets.py, remove_flagged_regions.py  -> $sample.in_good_regions.tsv
    bash helpers/remove_clusters.sh     $sample   # step 3  -> $sample.no_clustered_variants.tsv
    bash helpers/collapse_duplicates.sh $sample   # step 4  -> $sample.collapsed_duplicates.tsv
    bash helpers/remove_distortions.sh  $sample   # step 5  -> $sample.no_distortions.tsv
    # step 6: annotate_gnomad.py                   -> $sample.gnomad_annotated.tsv
    bash helpers/gnomad_filter.sh       $sample   # step 6  -> $sample.no_gnomad.tsv
    # steps 7-8: check_remapped_reads.py and single_tissue_combine_platform.py,
    #            once per tissue and platform     -> tissue_data/$sample.$tissue.tsv
done
bash helpers/combine_tissues.sh                   # step 9, all donors -> $sample.combined_tissues.tsv
for sample in SMHT001 SMHT004 SMHT005 SMHT012 SMHT024; do
    bash helpers/remove_tests.sh        $sample   # step 10 -> $sample.passed_strand_bias.tsv
    bash helpers/ga_filter.sh           $sample   # step 11 -> $sample.no_homopolymer.tsv (final call set)
done
```

`helpers/combine_tissues.sh` lists the tissues used for each donor. The final
per-donor call sets are the `$sample.no_homopolymer.tsv` files.

---

## License

Released under the MIT License (see `LICENSE`).

## Code availability

https://github.com/mdnoyes/smart_p25

# PyRelFit

`pyrelfit` is a command-line Python tool designed to compute the relative fitness of SNP (single-nucleotide polymorphism) mutations across multiple sample generations. It processes a VCF (Variant Call Format) file, splits variants by generation and chromosome, computes allele frequency ratios, and produces normalized relative fitness metrics. Optionally, it can generate scatter-plot graphics of normalized fitness values.

---

## Key Features

* **VCF Parsing**: Reads standard VCF files using `cyvcf2` with multi-threaded parsing support.
* **Generation-Based Splitting**: Dynamically groups samples into user-defined generations via regex patterns.
* **Chromosome Filtering**: Users can restrict analysis to a subset of chromosomes using regex.
* **Relative Fitness Computation**: Calculates allele frequency ratios and computes a fitness-like score `w` per SNP pair.
* **Normalization**: Scales the `RelativeFitness` values by the maximum observed per generation pair, constraining to \[0,1].
* **Graphics Support**: Optional scatter plots of normalized fitness values across chromosome positions.
* **Parallel Processing**: Leverages Python’s `multiprocessing.Pool` to parallelize VCF filtering, merging, and normalization across cores.
* **Optional Clean-up**: Temporary directories and files can be automatically removed or retained based on user preference.
* **Outlier Grouping**: Optionally collects high RelativeFitness variants into grouped CSV files by threshold (>0.8, >0.6, >0.4, >0.2) for further analysis.

---

## Dependencies

The following Python packages and tools are required:

* **Python 3.7+**
* [`cyvcf2`](https://github.com/brentp/cyvcf2): Fast VCF parsing.
* `click`: Command-line interface creation.
* `matplotlib`: Plotting library for generating graphics.
* Standard library modules: `csv`, `os`, `sys`, `time`, `re`, `glob`, `multiprocessing`.

You can install the Python dependencies via `pip`:

```bash
pip install cyvcf2 click matplotlib
```

Additional dependencies:
* `tabix`: Creation of `.tbi` index.
* `bgzip`: Compressing vcf file for faster processing.

*Both are part of `htslib`.*

---

## Installation

1. **Clone the repository**

2. **Install python dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

3. **Install `htslib`, required for indexing and compressing VCF files and `cyvcf`**

    * On Ubuntu/Debian:
    ```bash
    sudo apt-get update
    sudo apt-get install tabix
    ```

    * Via Conda (platform-independent):
    ```bash
    conda install -c bioconda htslib
    ```

    * On macOS (Homebrew):
    ```bash
    brew install htslib
    ```
---

## Usage

### Indexing VCF Files
cyvcf2 requires a Tabix index (.tbi) for random access to VCF/BCF files. You can create an indexed, compressed VCF using bgzip and tabix:

1. Compress the VCF
   ```bash
   bgzip -c input.vcf > input.vcf.gz
   ```

2. Create the VCF index
   ```bash
   tabix -p vcf input.vcf.gz
   ```
After indexing, ensure both input.vcf.gz and input.vcf.gz.tbi are in the working directory before running pyrelfit.

### Basic Command

```bash
python pyrelfit.py -i input.vcf -g /1/^gen1_/2/^gen2_/ -p 1_2 -c 4 -o results -t tmp
```

This command processes `input.vcf`, groups samples whose names match `^gen1_` into generation `1` and `^gen2_` into generation `2`, compares generation pair `1_2`, uses 4 CPU cores, writes outputs in `results/`, and uses `tmp/` for intermediate files.

### Options and Arguments

| Option                | Shortcut | Description                                                                                                                  | Default   |
|-----------------------|----------|------------------------------------------------------------------------------------------------------------------------------|-----------|
| `--input-file`        | `-i`     | Path to the input VCF file (required).                                                                                       |           |
| `--generations`       | `-g`     | Generation specification string: `/id1/regex1/id2/regex2/...` (required).                                                    |           |
| `--cores`             | `-c`     | Number of CPU cores to use for parallel processing.                                                                          | `1`       |
| `--chromosomes`       | `-C`     | Regex to filter chromosome names (e.g., `^chr[0-9]+$`).                                                                      | All       |
| `--generation-pairs`  | `-p`     | Comma-separated list of gen IDs to compare (e.g., `1_2,1_3`).                                                                | `1_3,2_3` |
| `--out-dir`           | `-o`     | Output directory for generated CSVs and optional plots.                                                                      | `results` |
| `--temp-dir`          | `-t`     | Temporary directory for intermediate CSVs.                                                                                   | `tmp`     |
| `--outliers`          | `-O`     | Flag to enable grouping and output of high relative fitness (RF) mutations into separate CSV files (>0.8, >0.6, >0.4, >0.2). | `False`   |
| `--generate-graphics` | `-G`     | Flag to enable generation of scatter-plot PNGs for each chromosome-pair.                                                     | `False`   |
| `--keep-temp`         | *None*   | Flag to retain temporary files after execution for debugging.                                                                | `False`   |

### Generation String Syntax

* Must begin and end with `/`, e.g., `/1/^genA_/2/^genB_/3/ctrl_/`.
* Each generation block consists of an **ID** (alphanumeric) followed by a **regex** to match sample names.
* IDs serve as keys in the output mapping and are used to name intermediate and result files.

### Examples

1. **Two generations (wildtype vs mutant)**

   ```bash
   python pyrelfit.py \
     --input-file experiment.vcf \
     --generations /wt/^WT_/mut/^MUT_/ \
     --generation-pairs wt_mut \
     --cores 8 \
     --generate-graphics true
   ```

   * Splits `WT_` samples into `wt`, `MUT_` into `mut`.
   * Compares relative fitness of mutants versus wildtype.
   * Uses 8 cores and generates PNG plots.

2. **Three generations (baseline, intermediate, final)**

   ```bash
   python pyrelfit.py -i lineage.vcf \
     -g /t0/^T0_/t1/^T1_/t2/^T2_/ \
     -p t0_t1,t1_t2 \
     -o lineage_results \
     -C "^chr[1-9]$" \
     -c 4
   ```

   * Only analyzes chromosome names `chr1`–`chr9`.

---

## Extending or Modifying

1. **Adding New Metrics**: The `compute_counts` and `process_pair` functions can be extended to compute alternative fitness metrics or summary statistics. Be cautious to adjust normalization accordingly.

2. **Alternative File Formats**: To support compressed VCFs (e.g., `.vcf.gz`), ensure the `VCF()` constructor and file open calls handle gzip files, or layer with Python’s `gzip` module.

3. **Custom Output**: Results are written as CSV per chromosome-pair. You can post-process these using R, Python, or Excel for downstream plotting or statistical tests.

---

## Internal Workflow

The pipeline proceeds in three major steps:

### Step 1: Filtering & Splitting

* **Function**: `filter_and_split`
* **Action**: Iterates through each specified chromosome, loads variants via `cyvcf2` limited to matching sample lists, and writes per-generation CSVs containing: position, REF, ALT, alt allele count, total depth.
* **Parallelization**: Each (chromosome, generation) pair is processed in parallel via `multiprocessing.Pool`.

### Step 2: Merging & Computation

* **Function**: `merge_and_compute`
* **Action**: For each generation pair (e.g. `1_2`), reads the two per-gen CSVs of each chromosome line-by-line, computes allele frequencies (`f1`, `f2`), then calculates a relative fitness weight `w = (f2^2) / (2*f1^2 - f1*f2^2)` under valid conditions.
* **Output**: Writes `chrom.<gen1>_<gen2>.csv` with columns `Pos`, `Ref`, `Alt`, `RF`.
* **Max Tracking**: Records maximum `RF` per generation pair for later normalization.

### Step 3: Normalization, Grouping & Plotting

* **Function**: `normalise`
* **Action**: Scales each `RF` value by the maximum observed for that pair, yielding values in `[0, 1]`. If graphics are enabled, generates scatter plots (`.png`) of normalized `RF` vs position.
* **Grouping**: When enabled, variants with high relative fitness (RF) values are categorized into groups and saved as separate CSV files. Group thresholds are: >0.8, >0.6, >0.4, and >0.2.
* **Graphics**: Uses `matplotlib` with parameters `GRAPHICS_OPACITY`, `GRAPHICS_POINT_SIZE`.

---

## Performance & Parallelization

* **IO Bound**: The pipeline reads/writes many small CSV files. Using an SSD and sufficient RAM cache improves throughput.
* **CPU Utilization**: VCF parsing (via `cyvcf2`) and CSV scanning dominate CPU load. The default single-core can be a bottleneck; set `--cores` to match available physical cores (avoid hyperthreading-only cores for best performance).
* **Memory Footprint**: Each VCF reader instantiation holds chromosome data per process. Avoid setting `--cores` too high on memory-limited systems.

---

## Troubleshooting & FAQs

**Q1: CSVs missing or empty?**

* Ensure sample regexes correctly match sample names in the VCF header. Run `cyvcf2.VCF(input.vcf).samples` interactively to inspect names.

**Q2: `ValueError: Invalid generation format`?**

* Generation string must start and end with `/`. It should contain an even number of segments: `/id1/regex1/id2/regex2/.../`.

**Q3: Graphics not generated?**

* Confirm you used `-G` or `--generate-graphics`. Ensure `matplotlib` is installed and you have write permissions to generate PNGs.

**Q4: Errors related to `htslib` linking?**

* Reinstall `cyvcf2` after installing system-level `htslib` headers (e.g., `sudo apt-get install libhts-dev`).

---

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.


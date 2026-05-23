# TEATIME
## esTimating EvolutionAry events Through sIngle-tiMepoint sEquencing

`TEATIME` analyses cancer sequencing samples based on variant allele frequency (VAF) data. It derives evolutionary parameters — mutation rate, selection coefficient, subclone emergence timing, and subclone expansion score — by examining the composition of VAF clusters.

📖 **[Example Workflow](https://htmlpreview.github.io/?https://github.com/liliulab/TEATIME/blob/main/vignettes/workflow.html)**

---

## Install

```r
# install.packages("devtools")
devtools::install_github("liliulab/TEATIME")
```

---

## Quick Start

TEATIME ships in two modes controlled by `fast_version`:

- **Fast mode (`fast_version = TRUE`)** — roughly **5–7× faster per sample** than default.
- **Default mode (`fast_version = FALSE`)** — the reference pipeline, may provide better estimation.

```r
library(TEATIME)
# input has 3 columns: REF, ALT, CN (per-mutation copy number)
result <- TEATIME.run(your_vcf_data, input_format = "vcf", beta = 0.9,
                      fast_version = TRUE)
print(result)
```

Cohort batch processing can layer `mclapply`  across samples on top of fast mode:

```r
parallel::mclapply(sample_files, function(f) {
  inp <- readRDS(f)
  TEATIME.run(inp, input_format = "magos", beta = 0.9,
              depth = round(mean(inp$result$depth.1)), fast_version = TRUE)
}, mc.cores = parallel::detectCores() - 1)
```

---

## Input Formats

TEATIME accepts three input formats via the `input_format` argument:

### Option 1 — Raw read counts (`input_format = "vcf"`)

A `data.frame` with three columns: reference counts, alternate counts, and copy number. Each row is one somatic mutation.

| REF | ALT | CN |
|-----|-----|-----|
| 120 |  30 |  2 |
| 200 |  50 |  2 |

```r
result <- TEATIME.run(your_vcf_data, input_format = "vcf", beta = 0.9)

# Optionally persist the intermediate MAGOS clustering for re-use later
TEATIME.run(your_vcf_data, input_format = "vcf", beta = 0.9,
            save_magos = TRUE)   # writes <output_folder>/<prefix>_MAGOS.rds
```

TEATIME runs MAGOS internally to cluster mutations before estimation. 


### Option 2 — MAGOS clustering result (`input_format = "magos"`)

A list with elements `purity` (numeric) and `result` (data.frame from MAGOS output):

```r
input <- list(purity = magos$purity, result = magos$results)
result <- TEATIME.run(input, input_format = "magos", beta = 0.9, depth = depth)
```

For more on MAGOS, see [github.com/liliulab/magos](https://github.com/liliulab/magos). An example of this procedure is provided in the [TEATIME example workflow](https://htmlpreview.github.io/?https://github.com/liliulab/TEATIME/blob/main/vignettes/workflow.html).

### Option 3 — Pre-clustered data (`input_format = "raw"`)

A `data.frame` with columns `vaf.1`, `depth.1`, and `colors` (cluster label):

```r
result <- TEATIME.run(my_data, input_format = "raw", beta = 0.9)
```

---

## Output

`TEATIME.run()` returns a one-row `data.frame`:

| Column | Description |
|--------|-------------|
| `name` | Sample ID |
| `mu` | Mutation rate |
| `s` | Selection coefficient |
| `emergence_time` | Emergence time of the subclone |
| `tau` | Subclone expansion score |
| `p` | Subclonal fraction |

When written to disk (`write_final = TRUE`), the file includes a `##` header line explaining each column.

---

## Key Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `input` | — | Input data (see Input Formats) |
| `beta` | `0.9` | Cell survival rate |
| `depth` | `NA` | Mean sequencing depth (auto-computed if `NA`) |
| `input_format` | `"vcf"` | One of `"vcf"`, `"magos"`, `"raw"` |
| `verbose` | `FALSE` | Print step-level progress |
| `output_folder` | `"./"` | Directory for output files |
| `output_prefix` | `"TEATIME"` | Prefix for output file names |
| `id` | `"T01"` | Sample identifier in result table |
| `write_final` | `TRUE` | Write `.final.txt` result file |
| `seed` | `123` | Random seed for reproducibility; set to `NA` to run estimators three times independently for stochastic robustness |
| `debug` | `FALSE` | Enable debug mode |
| `save_magos` | `FALSE` | `vcf` mode only: when `TRUE`, save the intermediate MAGOS clustering to `<output_folder>/<output_prefix>_MAGOS.rds`|
| `fast_version` | `FALSE` | When `TRUE`, vectorised for ~5–7× per-sample speed-up. |

---

## Debug Mode

Set `debug = TRUE` to trace exactly where the pipeline fails. Each step prints its name, key intermediate values, and elapsed time. On error, the failing step and error message are shown before stopping.

```r
TEATIME.run(input, debug = TRUE)
```

---

## Test Data

Two files are bundled in `inst/extdata` for validation:

| File | Description |
|------|-------------|
| `exampledata.rds` | Raw somatic mutation data (REF/ALT counts) |
| `MAGOS.rds` | Pre-computed MAGOS clustering result from the same data |

```r
magos <- readRDS(system.file("extdata", "MAGOS.rds", package = "TEATIME"))
input <- list(purity = magos$purity, result = magos$results)
result <- TEATIME.run(input, input_format = "magos", beta = 0.9, depth = round(mean(magos$results$depth.1)))
```

---

## Extensibility — Custom Growth Models 

TEATIME supports pluggable growth models. Register any model with `register_growth_model()`:

```r
# Model with extra parameters — declare ctx to receive the pipeline context
register_growth_model("logistic", function(i, p, beta, ctx) {
  K <- ctx$extra$carrying_capacity  # any extra param passed via extra=list(...)
  p / 2 + (1 - p) / (2 * (1 + (K / p - 1) * exp(-beta * i)))
})

result <- TEATIME.run(input, growth_model = "logistic",
                      extra = list(carrying_capacity = 2))
```

The function must accept `i`, `p`, and `beta` as named arguments and return a single numeric VAF value. Declare a `ctx` argument if your model needs additional parameters — TEATIME will pass the full pipeline context, and you can store extras in `ctx$extra` via the `extra` argument to `TEATIME.run()`.

---

## Reference

Chen H, Shu J, Mudappathi R, Li E, Wang P, Bergsagel L, Yang P, Sun Z, Zhao L, Shi C, Townsend JP, Maley C, Liu L. Competing subclones and fitness diversity shape tumor evolution across cancer types. Bioinformatics. 2026 Feb 28;42(3):btag127. doi: 10.1093/bioinformatics/btag127. PMID: 41826799; PMCID: PMC13025073.

## Contributors

Algorithm of TEATIME was developed by Li Liu and Hai Chen. Please contact liliu at asu.edu for questions or suggestions.

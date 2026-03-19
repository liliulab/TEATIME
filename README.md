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

**Dependencies**: `dplyr`, `igraph`, `Matrix`, `strucchange`, `likelihoodExplore`, `RBesT (>= 1.6.6)`, `magrittr`

---

## Quick Start

```r
library(TEATIME)

# Load the bundled MAGOS clustering result
magos <- readRDS(system.file("extdata", "MAGOS.rds", package = "TEATIME"))
input <- list(purity = magos$purity, result = magos$results)
depth <- round(mean(magos$results$depth.1))

result <- TEATIME.run(input, beta = 0.9, depth = depth, seed = 123)
print(result)
#   name   mu      s  emergence_time   tau       p
#   T01   17.9  0.442       10        3.12   0.406
```

---

## Input Formats

TEATIME accepts three input formats via the `input_format` argument:

### Option 1 — Raw read counts (`input_format = "vcf"`)

A `data.frame` with three columns: reference counts, alternate counts, and copy number. Each row is one somatic mutation.

| REF | ALT | CN |
|-----|-----|----|
| 120 |  30 |  2 |
| 200 |  50 |  2 |

```r
# e.g. your_vcf_data is a data.frame with REF, ALT, CN columns
result <- TEATIME.run(your_vcf_data, input_format = "vcf", beta = 0.9, seed = 123)
```

TEATIME runs MAGOS internally to cluster mutations before estimation.

### Option 2 — MAGOS clustering result (`input_format = "magos"`, default)

A list with elements `purity` (numeric) and `result` (data.frame from MAGOS output):

```r
input <- list(purity = magos.33$purity, result = magos.33$results)
result <- TEATIME.run(input, beta = 0.9, depth = 1000, seed = 123)
```

For more on MAGOS, see [github.com/liliulab/magos](https://github.com/liliulab/magos).

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
| `mu` | Mutation rate (mutations per cell division) |
| `s` | Selection coefficient |
| `emergence_time` | Emergence time of the subclone (cell divisions) |
| `tau` | Subclone expansion score (`tend / emergence_time`) |
| `p` | Clonal fraction |

When written to disk (`write_final = TRUE`), the file includes a `##` header line explaining each column.

---

## Key Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `input` | — | Input data (see Input Formats) |
| `beta` | `0.9` | Cell survival rate |
| `depth` | `NA` | Mean sequencing depth (auto-computed if `NA`) |
| `p_thre` | `0.01` | P-value threshold for breakpoint tests |
| `input_format` | `"magos"` | One of `"magos"`, `"vcf"`, `"raw"` |
| `growth_model` | `"exponential"` | Growth model name (see Extensibility) |
| `verbose` | `FALSE` | Print step-level progress |
| `output_folder` | `"./"` | Directory for output files |
| `output_prefix` | `"TEATIME"` | Prefix for output file names |
| `id` | `"T01"` | Sample identifier in result table |
| `write_final` | `TRUE` | Write `.final.txt` result file |
| `seed` | `NA` | Random seed for reproducibility |

---

## Test Data

Three files are bundled in `inst/extdata` for validation:

| File | Description |
|------|-------------|
| `exampledata.rds` | Raw somatic mutation data (REF/ALT counts) |
| `MAGOS.rds` | Pre-computed MAGOS clustering result from the same data |
| `TEATIME.final.txt` | Reference output from the original TEATIME pipeline |

```r
# Access with system.file()
magos <- readRDS(system.file("extdata", "MAGOS.rds", package = "TEATIME"))
input <- list(purity = magos$purity, result = magos$results)
result <- TEATIME.run(input, beta = 0.9,
                      depth = round(mean(magos$results$depth.1)),
                      seed = 123)
```

---

## Extensibility — Custom Growth Models

TEATIME v2 supports pluggable growth models. Register any model with `register_growth_model()`:

```r
# Simple model — just needs i, p, beta
register_growth_model("exponential2", function(i, p, beta) {
  p / 2 + (1 - p) / (2 * exp(log(2) * beta * i))
})

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

Manuscript in preparation.

## Contributors

Algorithm of TEATIME was developed by Li Liu and Hai Chen. Please contact liliu at asu.edu for questions or suggestions.

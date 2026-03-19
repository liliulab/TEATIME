# TEATIME Example Workflow

*2025-10-11*

---

## Overview

TEATIME is a framework for estimating evolutionary parameters using single-timepoint WES or WGS data. In this workflow, we demonstrate the complete process using an example from the MOBSTER dataset, including data preprocessing, clustering, and TEATIME analysis.

---

## Example Case

This example contains a simulated single-timepoint tumor composed of a root clone (`clone0`) and one expanding subclone (`clone1`). It includes the overall cellular composition, evolutionary parameters, and simulated sequencing reads.

---

## Clonal Composition

```r
library(dplyr)

simulations <- readRDS(system.file("extdata", "exampledata.rds", package = "TEATIME"))

counts <- simulations$cell_counts["counts", ]
clone0_fraction <- as.numeric(counts["clone0"] / counts["total"])
clone1_fraction <- as.numeric(counts["clone1"] / counts["total"])

cat("Clone0 fraction (within tumor):", round(clone0_fraction * 100, 1), "%\n")
## Clone0 fraction (within tumor): 55.7 %
cat("Clone1 fraction (within tumor):", round(clone1_fraction * 100, 1), "%\n")
## Clone1 fraction (within tumor): 44.3 %
```

The two clones are roughly balanced (~56% vs ~44%) in this simulation.

---

## Evolutionary Parameters

We examine the true evolutionary parameters — mutation rate (μ) and relative fitness (s) — used in the simulation.

```r
evo.parameter <- simulations$clone_parameters
mu <- evo.parameter$mutation_rate[1]
s  <- (evo.parameter$birth_rate[2] - evo.parameter$death_rates[2]) /
      (evo.parameter$birth_rate[1] - evo.parameter$death_rates[1]) - 1

cat("Mutation rate:", mu, "\n")
## Mutation rate: 16
cat("s (clone1 vs clone0):", round(s, 2), "\n")
## s (clone1 vs clone0): 0.5
```

---

## Sequencing Data

We examine the sequencing data providing VAFs for parameter inference. Since this simulation assumes a diploid genome and 100% purity, no additional filtering is required.

> **Note:** In real datasets, retain only mutations in copy-neutral (CN = 2) regions, as TEATIME assumes diploid copy number.

```r
library(data.table)

data <- as.data.table(simulations$sequencing)
head(data)
##       VAF   ALT    DP CLONE TRUE_VAF TRUE_CLUSTER
## 1: 0.0676     5    74     1   0.0278           19
## 2: 0.0526     6   114     1   0.0226           1c
## 3: 0.0530     7   132     1   0.0568           15
## 4: 0.0609     7   115     1   0.0568           15
## 5: 0.0568     5    88     1   0.0568           15
## 6: 0.0588     6   102     1   0.0568           15

data[, CN := 2]
```

---

## Analysis

### Input Options

TEATIME performs inference based on clustering results. Three options are available:

1. **Raw read counts** — TEATIME runs MAGOS internally (`input_format = "vcf"`)
2. **MAGOS result** — run MAGOS externally and pass its output (`input_format = "magos"`)
3. **Custom clustering** — any method (e.g. PyClone), formatted as `vaf.1 / depth.1 / colors` (`input_format = "raw"`)

---

### Option 2 — MAGOS Clustering (Recommended)

Prepare REF/ALT input and run MAGOS. Save the result for reproducibility:

```r
library(MAGOS)

data[, t_ref_count := DP - ALT]
data[, t_alt_count := ALT]
input.data <- data[, .(REF = t_ref_count, ALT = t_alt_count)]

# Run MAGOS and save for reproducibility
run <- mag.single.run(input.data = input.data, fold = TRUE)
saveRDS(run, "MAGOS.rds")
```

Load the pre-computed result and inspect clusters:

```r
input.file <- readRDS(system.file("extdata", "MAGOS.rds", package = "TEATIME"))

plot(input.file$results$vaf.1, input.file$results$depth.1,
     col  = input.file$results$colors,
     pch  = 19, cex = 0.6,
     xlab = "VAF", ylab = "Depth",
     main = "MAGOS Clustering Result")
```

MAGOS identified three clusters, indicating the presence of a subclone. This provides TEATIME with the necessary input to estimate evolutionary parameters.

---

## Run TEATIME

Pass the MAGOS result directly to `TEATIME.run()`. Sequencing depth is optional — TEATIME computes it from the data if omitted.

```r
library(TEATIME)

samplename    <- "test"
output.folder <- file.path("./", samplename)
dir.create(output.folder, showWarnings = FALSE)

depth <- round(mean(input.file$results$depth.1))  # optional

TEATIME.result <- TEATIME.run(
  input         = input.file,
  input_format  = "magos",
  beta          = 0.9,
  depth         = depth,
  output_folder = output.folder,
  output_prefix = "TEATIME",
  id            = samplename,
  seed          = 123
)

print(TEATIME.result)
##   name   mu        s  emergence_time      tau         p
## 1 test 17.9 0.440502            10  3.131972  0.4062049
```

### Output Columns

| Column | Description |
|--------|-------------|
| `name` | Sample ID |
| `mu` | Mutation rate (mutations per cell division) |
| `s` | Selection coefficient — relative growth advantage of the subclone |
| `emergence_time` | Emergence time of the subclone (cell divisions) |
| `tau` | Subclone expansion score (tend / emergence_time) |
| `p` | Cellular frequency of the subclone |

The output file `TEATIME.final.txt` includes a `##` header line explaining each column.

### Compare to True Values

```r
mu_error <- (TEATIME.result$mu - mu) / mu
s_error  <- (TEATIME.result$s  - s)  / s
p_error  <- (TEATIME.result$p  - clone1_fraction) / clone1_fraction

cat(sprintf("Relative error of mutation rate: %.2f%%\n", mu_error * 100))
## Relative error of mutation rate: 11.88%
cat(sprintf("Relative error of fitness:       %.2f%%\n", s_error  * 100))
## Relative error of fitness:       -11.90%
cat(sprintf("Relative error of frequency:     %.2f%%\n", p_error  * 100))
## Relative error of frequency:     -8.30%
```

---

## Alternative Input Options

### Option 1 — Raw Sequencing Data

Pass the read count table directly; TEATIME runs MAGOS internally.

```r
input.data <- data[, .(REF = DP - ALT, ALT = ALT, CN = CN)]

TEATIME.run(
  input         = input.data,
  input_format  = "vcf",
  output_folder = output.folder,
  id            = samplename
)
```

### Option 3 — Custom Clustering

Format any clustering result as a `data.frame` with columns `vaf.1`, `depth.1`, and `colors`. Column names must be exact.

```r
input.data <- data.frame(
  vaf.1   = data$VAF,
  depth.1 = data$DP,
  colors  = your_cluster_labels   # from PyClone or any other tool
)

TEATIME.run(
  input         = input.data,
  input_format  = "raw",
  output_folder = output.folder,
  id            = samplename
)
```

> **Fallback:** When diploid-region mutations are too few, tools such as PureCN can provide adjusted VAFs. Simulate REF/ALT counts from the mean depth and adjusted VAF, then use `input_format = "vcf"`. This is a last resort — interpret results with caution.

---

## Custom Growth Models

TEATIME v2 supports pluggable growth models via `register_growth_model()`. Declare an optional `ctx` argument in your function to receive extra parameters passed through `extra = list(...)`.

### Logistic Growth

Accounts for a finite carrying capacity — useful when resource competition limits tumour expansion:

```r
register_growth_model("logistic", function(i, p, beta, ctx) {
  K <- ctx$extra$carrying_capacity %||% 1   # carrying capacity, default 1
  p / 2 + (1 - p) / (2 * (1 + (K / p - 1) * exp(-beta * i)))
})

TEATIME.run(
  input         = input.file,
  input_format  = "magos",
  growth_model  = "logistic",
  extra         = list(carrying_capacity = 2),
  output_folder = output.folder,
  id            = samplename
)
```

---

## References

Caravagna, G., Heide, T., Williams, M.J. et al. Subclonal reconstruction of tumors by using machine learning and population genetics. *Nat Genet* 52, 898–907 (2020). <https://doi.org/10.1038/s41588-020-0675-5>

Navid Ahmadinejad et al. Accurate Identification of Subclones in Tumor Genomes, *Molecular Biology and Evolution*, Volume 39, Issue 7, July 2022, msac136. <https://doi.org/10.1093/molbev/msac136>

Riester, M. et al. PureCN: copy number calling and SNV classification using targeted short read sequencing. *Source Code Biol Med* 11, 13 (2016). <https://doi.org/10.1186/s13029-016-0060-z>

Roth, A. et al. PyClone: statistical inference of clonal population structure in cancer. *Nat Methods* 11, 396–398 (2014). <https://doi.org/10.1038/nmeth.2883>

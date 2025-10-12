# TEATIME
## esTimating EvolutionAry events Through sIngle-tiMepoint sEquencing

 `TEATIME` analyzes cancer sequencing samples based on allele frequency data. It derives evolutionary parameters by examining the composition of clusters formed by various VAFs (Variant Allele Frequencies). These parameters include mutation rate, subclone fitness, the timing of driver mutation occurrence, and the timing of sample collection.
 

## Install

```
library(devtools)
devtools::install_github('liliulab/TEATIME')
```


## Input Data
### TEATIME provides three options for input files/data.

Option 1: Reference & Alternate Read Counts with CN 

`TEATIME` can utilize reference and alternate read counts directly from a VCF file. In the standard VCF format, this information is found in the 'AD' (Allelic Depths) section of the FORMAT column.

The input data for single sample analysis should be a dataframe or a matrix with three columns. Each row corresponds to a point mutation. The reference reads in the first column and the alterante reads in the second column, Copy number information in third column.

| | Ref Counts | Alt Counts| CN |
|----|---------|--------|--------|
|mut 1  | ref 1      | alt 1     | 2 |
|mut 2 | ref 2      | alt 2     | 2 |
|mut 3 | ref 3      | alt 3     | 2 |
|...|...|...|


To run `TEATIME` using this input, include `0` in `steps`. `TEATIME` will run `MAGOS` in this case.

#### Usage
```
TEATIME.run(input.file,steps=0:5)
```
<br/><br/>


Option 2: MAGOS result

The second option is to run `MAGOS` on your local side, then use the returned result from `MAGOS` as the input for `TEATIME`. In this case, set `magos_object = TRUE` and ensure `0` is not included in the `steps`.

#### Usage
```
TEATIME.run(input.file,magos_object = TRUE,steps=1:5)
```
For more information regarding `MAGOS`, visit the [MAGOS](https://github.com/liliulab/magos/).
<br/><br/>


Option 3: Own Sequencing Data Clustering Result

The third option is to run TEATIME with your own sequencing data clustering result. The input must be a dataframe with three columns (`vaf.1, depth.1, colors`). Each row corresponds to a point mutation.The first two columns are straightforward, while the third column, `colors`, refers to the cluster corresponding to each mutation. Do not change the column names.

Note that in this case, set `magos_object = FALSE` and ensure `0` is not included in the `steps`.

| vaf.1 | depth.1 | colors |
|----|---------|--------|
|vaf 1  | dp1     | c1    |
|vaf 2 | dp2    | c2     |
|vaf 3 | dp3   | c3   |
|...|...|...|

#### Usage
```
TEATIME.run(input.file,magos_object = FALSE,steps=1:5)
```
<br/><br/>

An example of this procedure is provided in the [TEATIME example workflow](https://haichen294.github.io/teatime-workflow/).

Additional test data can be found in the data folder. To ensure the reproducibility of the test and control the consistency of the simulation, set `seed=123` (or any number). `123` is the default seed for generating `magos.sample.rda`. 
This setting is limited to repeating the same results and is not recommended for real analysis.

#### Test code
```
TEATIME.run(input.file=vcf.33,depth=1000,steps=0:5,seed = 123) 

TEATIME.run(input.file=magos.33,depth=1000,magos_object=T,steps=1:5,seed = 123) 
```



## Reference

Please see preprint: https://www.biorxiv.org/content/10.1101/2025.05.31.657191v1
<br/><br/>


## Contributors
Algorithm of TEATIME was developed by Li Liu and Hai Chen. Please contact liliu at asu.edu for any questions or suggestions. 

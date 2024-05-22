# TEATIME
## esTimating EvolutionAry events Through sIngle-tiMepoint sEquencing

 TEATIME analyzes cancer sequencing samples based on allele frequency data. It derives evolutionary parameters by examining the composition of clusters formed by various VAFs (Variant Allele Frequencies). These parameters include mutation rate, subclone fitness, the timing of driver mutation occurrence, and the timing of sample collection.
 

## Install

```
library(devtools)
devtools::install_github('liliulab/MAGOS')
```


## Input Data
MAGOS only requires the reference and the alternate read counts from the VCF file. In the the standard VCF format this information can be extracted from the 'AD' (Allelic Depths) section in the FORMAT column. 

For more information, visit the [MAGOS homepage](https://github.com/liliulab/magos/).

### Single Sample
The input data for single sample analysis should be a dataframe or a matrix with two columns. Each row corresponds to a point mutation. The reference reads in the first column and the alterante reads in the second column. 

| | Ref Counts | Alt Counts|
|----|---------|--------|
|mut 1  | ref 1      | alt 1     |
|mut 2 | ref 2      | alt 2     |
|mut 3 | ref 3      | alt 3     |
|...|...|...|


| | Ref Sample 1 | Alt Sample 1| Ref Sample 2| Alt Sample 2| Ref Sample 3 | Alt Sample 3|...|...|
|----|---------|--------|---|---|---|---|---|---|
|mut 1  | ref 1 s1     | alt 1 s1    | ref 1 s2| alt 1 s2| ref 1 s3 | alt 1 s3|...|...|
|mut 2 | ref 2 s1     | alt 2 s1     |ref 2 s2| alt 2 s2| ref 2 s3 | alt 2 s3|...|...|
|mut 3 | ref 3 s1    | alt 3 s1    |ref 3 s2| alt 3 s2| ref 3 s3 | alt 3 s3|...|...|
|...|...|...|...|...|...|...|...|...|


## Usage
```
run = mag.single.run(input.data)


# the results will be in run$results
# the results can be visualized using the following script: 

plot(run$results$vaf.1, run$results$depth.1, col= run$results$colors, xlab= 'VAF', ylab='Depth') 
```

Examples and test data and explanations on different functionalities can be found in the ACE_workshop folder. 




## Reference

The paper in under manuscript. 




## Contributors
Algorithm of MAGOS was developed by Li Liu and Hai Chen. Please contact liliu at asu.edu for any questions or suggestions. 

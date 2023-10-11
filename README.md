# INCAData

## Installation
If the installation fails, install the dependencies:
```R
## data.table
install.packages("data.table")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## GenomicRanges
BiocManager::install("GenomicRanges")

## Rsamtools
BiocManager::install("Rsamtools")
```

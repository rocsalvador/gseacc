# gseacc

## R 

### Package installation

#### Requirements

- Rcpp (R package)
- g++

```bash
git clone https://github.com/rocsalvador/gseacc/
cd gseacc
make install -j
```

### Documentation

For the R documentation use ? in R. Available commands:

- ```?Gsea```
- ```?run```
- ```?runChunked```
- ```?filterResults```
- ```?normalizeExprMatrix```
- ```?readCsv```
- ```?readGeneSets```
- ```?writeGeneSets```

Example R scripts in [Efficient-rank-based-statistic-for-partially-overlapping-genesets](https://github.com/rocsalvador/Efficient-rank-based-statistic-for-partially-overlapping-genesets)

## C++

### How to run

Clone and compile:

```bash
git clone https://github.com/rocsalvador/gseacc/
cd gseacc
make cc
```
Create a config file `gsea.config` in the same folder where you will run the executable:

```bash
expression-matrix-file:     relative or full path to expression matrix csv file
sep:                        csv element separator for the expression matrix csv (t for tabular)
gene-sets-file:             relative or full path to the gene sets file
sep:                        csv element separator for the gene sets file (t for tabular)
output-file:                relative or full path to the GSEA results csv
sep:                        csv element separator for the GSEA results csv (t for tabular)
threads-used:               threads used for the GSEA computation (0 to use all available threads)
normalized-data:            0 if data is not normalized, 1 if it is not
ioutput:                    number of gene sets between std output
scrna:                      0 if it is a rna experiment (runRna), 1 if it is a sc-rna experiment (runScRna)
batch-size:                 number of lines read every loop for runScRna function 
```

Run the executable (`gsea.config` must be in the same folder from where you run the executable):

```bash
./gsea
```

### Documentation

For the C++ functions documentation check this [page](https://rocsalvador.github.io/).

## About

This package is part of a bachelor's thesis written at NTNU and supervised by Pål Sætrom (pal.satrom@ntnu.no).

### Author

Roc Salvador Andreazini (roc.salvador@estudiantat.upc.edu) 


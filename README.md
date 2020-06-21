# MIS-C Manuscript analysis and figures

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)
* [Omics integration](#omics-integration)
* [Setup](#setup)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink - NPX values)
- FACS
- Autoantibodies
	
## Dependencies
Project is created with:
* RStudio version: 3.6.2

## Repo description
- ```Figure*``` reproducible code for figures laid out in manuscript

## Omics integration
### MOFA
- Multi-Omics Factor Analysis (MOFA+) was used in this study in order to deconvolute the main sources of variation in the differents sets of data mentioned above. For more information, read their published Methods paper [Argelaguet et al. (2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1). 
- MOFA is publicly accessible here: https://github.com/bioFAM/MOFA2

## Setup
To run this project, install it locally using devtools:

```
$ install.packages('devtools')
$ library(devtools)
$ install_github('Brodinlab/MIS-C_manuscript')
$ library(MIS-C_manuscript)
```

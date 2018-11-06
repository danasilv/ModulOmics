# ModulOmics

## Summary
ModulOmics is a method to identify cancer driver pathways de novo by integrating multiple data types (protein-protein interactions, mutual exclusivity of mutations or copy number alterations, transcriptional co-regulation, and RNA co-expression) into a single probabilistic model. ModulOmics uses patient data, and protein-protein interaction networks, based here on the Hippie database and regulatory connections, based here on the Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining database (TRRUST).

## Installation
ModulOmics is available as ```R``` code. It has been tested on Mac and Linux, on R 3.4.3. It needs an active ```cplex``` (IBM CPLEX Optimizer) instance installed, and depends on the following R packages:
* ```TiMEx```, current version downloadable from [its github repository](https://github.com/csimona/TiMEx/)
* ```igraph``` v. 1.1.2
* ```cplexAPI``` v. 1.3.3
* ```gtools``` v. 3.5.0

## Running via snakefile
To install and run via snakefile please download the directory to your computer, and run using "snakemake".
The software requires an active ```cplex``` (IBM CPLEX Optimizer) instance installed, and the snakefile requires requires Python 3.5+ with snakemake installed (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). This will both install dependencies and the softare and run ModulOmics with the default parameters.

## Static data
The static data used in the MopdulOmics paper is available in the webserver link: (http://anat.cs.tau.ac.il/ModulOmicsServer/) and include:
1. Regulatory connections, based here on the TRRUST database (Han et. al. Scientific Reports, 2015).
2. Shortest paths as computed on top of the Hippie protein-protein interaction network (Schaefer et. al. PLoS ONE, 2012).
The dynamic portion of the data used in the ModulOmics paper, namely the genetic alterations and the gene expression per patient, is available in TCGA (https://cancergenome.nih.gov/).

For more details on the input and output structures of ModulOmics, please see the pdf manual.


